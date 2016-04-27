"""
This module contains functions which create NIRSPEC WCS reference files
in ASDF format from reference files delivered by the NIRSPEC team in
various other formats.

Note:
The team delivered reference files call forward transform the transform
from sky to detector. In the WCS pipeline this is the backward transform.
All ASDF reference files take this into account.

CDP4 delivery
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import glob
import os.path
import numpy as np
import re
import string
from astropy.io import fits
from asdf import AsdfFile
from astropy.modeling import models
from astropy.modeling.models import Mapping, Identity


def common_reference_file_keywords(reftype, title):
    """Create common reference file keywords."""
    ref_file_common_keywords = {"reftype": reftype,
                                "title": title,
                                "pedigree": "GROUND",
                                "author": "N. Dencheva",
                                "instrument": "NIRSPEC",
                                "exp_type": "NRS_FIXEDSLIT"
                                }
    return ref_file_common_keywords


def dict_from_file(flist):
    """Read in a file and return a dict of the key value pairs.

    This is a generic read for a text file with the following format

    keywords which start with any punctuation character except + or -
    keyword sections separated by a blank line

    example file: disperser_PRISM_TiltX.gtp
    if you feed a list of files the last unique parameter from all the files
    will be saved. Multiple values are saved as a list
    """
    # + and - not allowed as tokens
    tokens = string.punctuation.replace('-', '').replace('+', '')
    letters = re.compile("(^[a-zA-Z])")  # starts with a letter, make string
    numbers = re.compile("(^(?:[+\-])?(?:\d*)(?:\.\d*)?(?:[eE][+\-]?\d*$)?)")
    empty = re.compile("(^\s*$)")  # is a blank line

    if isinstance(flist, str):
        flist = [flist]
    parameters = dict()
    for name in flist:
        print("\nReading {}  ...".format(name))
        with open(name, 'r') as fh:
            lines = fh.readlines()
        value = list()
        key = ''

        for line in lines:
            if not empty.match(line):
                if line[0] in tokens:
                    if not key:
                        key = line[1:].split()[0].strip()  # first word only
                elif letters.match(line):
                    value.append(line.strip())
                elif numbers.fullmatch(line):
                        value.append(eval(line))
                elif numbers.match(line):
                    if numbers.fullmatch(line):
                        value.append(eval(line.strip()))
                    else:
                        newline = line.split()
                        if len(newline) > 0:
                            for n in newline:
                                if numbers.fullmatch(n):
                                    value.append(eval(n))
                                else:
                                    value.append(n.strip())
                else:
                    value.append(line.strip())
            else:
                if key and value:
                    parameters[key] = value
                    print("Setting {0:s} = {1}".format(key, value))
                key = ''
                value = list()
        if key and value:  # save final value if file ends without blank lines
            parameters[key] = value
            print("Setting {0:s} = {1}".format(key, value))

        return parameters


def coeffs_from_pcf(degree, coeffslist):
    coeffs = {}
    k = 0
    for i in range(degree + 1):
        for j in range(degree + 1):
            if i+j < degree+1:
                name = "c{0}_{1}".format(i, j)
                coeffs[name] =  float(coeffslist[k])
                k +=1
            else:
                continue
    return coeffs




def pcf_forward(pcffile, outname):
    """
    Create the **IDT** forward transform from collimator to gwa.
    """
    with open(pcffile) as f:
        lines = [l.strip() for l in f.readlines()]

    factors = lines[lines.index('*Factor 2') + 1].split()
    # factor==1/factor in backward msa2ote direction and factor==factor in sky2detector direction
    scale = models.Scale(float(factors[0]), name="x_scale") & \
          models.Scale(float(factors[1]), name="y_scale")

    rotation_angle = lines[lines.index('*Rotation') + 1]
    # The minius sign here is because astropy.modeling has the opposite direction of rotation than the idl implementation
    rotation = models.Rotation2D(-float(rotation_angle), name='rotation')


    # Here the model is called "output_shift" but in the team version it is the "input_shift".
    input_rot_center = lines[lines.index('*InputRotationCentre 2') + 1].split()
    input_rot_shift = models.Shift(-float(input_rot_center[0]), name='input_x_shift') & \
                 models.Shift(-float(input_rot_center[1]), name='input_y_shift')


    # Here the model is called "input_shift" but in the team version it is the "output_shift".
    output_rot_center = lines[lines.index('*OutputRotationCentre 2') + 1].split()
    output_rot_shift = models.Shift(float(output_rot_center[0]), name='output_x_shift') & \
                  models.Shift(float(output_rot_center[1]), name='output_y_shift')

    degree = int(lines[lines.index('*FitOrder') + 1])
    xcoeff_index = lines.index('*xForwardCoefficients 21 2')
    xlines = lines[xcoeff_index + 1: xcoeff_index + 22]
    xcoeff_forward = coeffs_from_pcf(degree, xlines)
    x_poly_forward = models.Polynomial2D(degree, name='x_poly_forward', **xcoeff_forward)

    ycoeff_index = lines.index('*yForwardCoefficients 21 2')
    ycoeff_forward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
    y_poly_forward = models.Polynomial2D(degree, name='y_poly_forward', **ycoeff_forward)

    xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
    xcoeff_backward = coeffs_from_pcf(degree, lines[xcoeff_index + 1: xcoeff_index + 22])
    x_poly_backward = models.Polynomial2D(degree, name='x_poly_backward', **xcoeff_backward)

    ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
    ycoeff_backward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
    y_poly_backward = models.Polynomial2D(degree, name='y_poly_backward', **ycoeff_backward)

    x_poly_forward.inverse = x_poly_backward
    y_poly_forward.inverse = y_poly_backward

    poly_mapping1  = Mapping((0, 1, 0, 1))
    poly_mapping1.inverse = Identity(2)
    poly_mapping2 = Identity(2)
    poly_mapping2.inverse = Mapping((0, 1, 0, 1))

    model = input_rot_shift | rotation | scale | output_rot_shift | \
          poly_mapping1 | x_poly_forward & y_poly_forward | poly_mapping2
    f = AsdfFile()
    f.tree = {'model': model}
    f.write_to(outname)


def pcf2asdf(pcffile, outname, ref_file_kw):
    """
    Create an asdf reference file with the transformation coded in a NIRSPEC
    Camera.pcf or Collimator*.pcf file.

    - forward (team): sky to detector
      - Shift inputs to input_rotation_center
      - Rotate inputs
      - Scale  inputs
      - Shift inputs to output_rot_center
      - Apply polynomial distortion
    - backward_team (team definition) detector to sky
      - Apply polynomial distortion
      - Shift inputs to output_rot_center
      - Scale  inputs
      - Rotate inputs
      - Shift inputs to input_rotation_center

    WCS implementation
    - forward: detector to sky
      - equivalent to backward_team
    - backward: sky to detector
      - equivalent to forward_team

    Parameters
    ----------
    pcffile : str
        one of the NIRSPEC ".pcf" reference files provided by the IDT team.
        "pcf" stands for "polynomial coefficients fit"
    outname : str
        Name of reference file to be wriiten to disk.

    Returns
    -------
    fasdf : AsdfFile
        AsdfFile object

    Examples
    --------
    >>> pcf2asdf("Camera.pcf", "camera.asdf")

    """
    linear_det2sky = linear_from_pcf_det2sky(pcffile)

    with open(pcffile) as f:
        lines = [l.strip() for l in f.readlines()]

    degree = int(lines[lines.index('*FitOrder') + 1])

    xcoeff_index = lines.index('*xForwardCoefficients 21 2')
    xlines = lines[xcoeff_index + 1: xcoeff_index + 22]
    xcoeff_forward = coeffs_from_pcf(degree, xlines)
    x_poly_backward = models.Polynomial2D(degree, name='x_poly_backward', **xcoeff_forward)

    ycoeff_index = lines.index('*yForwardCoefficients 21 2')
    ycoeff_forward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
    y_poly_backward = models.Polynomial2D(degree, name='y_poly_backward', **ycoeff_forward)

    xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
    xcoeff_backward = coeffs_from_pcf(degree, lines[xcoeff_index + 1: xcoeff_index + 22])
    x_poly_forward = models.Polynomial2D(degree, name='x_poly_forward', **xcoeff_backward)

    ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
    ycoeff_backward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
    y_poly_forward = models.Polynomial2D(degree, name='y_poly_forward', **ycoeff_backward)

    x_poly_forward.inverse = x_poly_backward
    y_poly_forward.inverse = y_poly_backward

    output2poly_mapping = Identity(2, name='output_mapping')
    output2poly_mapping.inverse = Mapping([0, 1, 0, 1])
    input2poly_mapping = Mapping([0, 1, 0, 1], name='input_mapping')
    input2poly_mapping.inverse = Identity(2)

    model_poly = input2poly_mapping | (x_poly_forward & y_poly_forward) | output2poly_mapping

    model = model_poly | linear_det2sky

    f = AsdfFile()
    f.tree = ref_file_kw.copy()
    f.tree['model'] = model
    f.write_to(outname)


def homothetic_sky2det(input_center, angle, scale, output_center):
    """
    Create the homothetic transform from a .pcf file.

    The forward direction is sky to detector.
    Parameters
    ----------
    input_center : ndarray of shape (1,2) or iterable
        (x, y) coordinate of the input rotation center
    angle : float
        Rotation angle in degrees
    scale : ndarray of shape (1,2) or iterable
        Scaling factors in x, y directions
    output_center : ndarray of shape (1,2) or iterable
        (x, y) coordinate of the output rotation center

    """

    input_center = np.array(input_center, dtype=np.float)
    output_center = np.array(output_center, dtype=np.float)
    scale = np.array(scale, dtype=np.float)
    angle = np.deg2rad(angle)

    rotmat_sky2det = [[np.cos(angle), np.sin(angle)],
                      [-np.sin(angle), np.cos(angle)]]
    scaling = np.array([[scale[0] , 0], [0, scale[1]]])
    mat_sky2det = np.dot(scaling, rotmat_sky2det)

    aff = models.AffineTransformation2D(matrix=mat_sky2det)


    transform = models.Shift(-input_center[0]) & models.Shift(-input_center[1]) | aff | \
              models.Shift(output_center[0]) & models.Shift(output_center[1])
    return transform


def homothetic_det2sky(input_center, angle, scale, output_center):
    """
    Create the homothetic transform from a .pcf file.

    The forward direction is sky to detector.
    Parameters
    ----------
    input_center : ndarray of shape (1,2) or iterable
        (x, y) coordinate of the input rotation center
    angle : float
        Rotation angle in degrees
    scale : ndarray of shape (1,2) or iterable
        Scaling factors in x, y directions
    output_center : ndarray of shape (1,2) or iterable
        (x, y) coordinate of the output rotation center

    """

    input_center = np.array(input_center, dtype=np.float)
    output_center = np.array(output_center, dtype=np.float)
    scale = np.array(scale, dtype=np.float)
    angle = np.deg2rad(angle)

    rotmat_sky2det = [[np.cos(angle), -np.sin(angle)],
                      [np.sin(angle), np.cos(angle)]]
    scaling = np.array([[1/scale[0] , 0], [0, 1/scale[1]]])
    mat_sky2det = np.dot(rotmat_sky2det, scaling)

    aff = models.AffineTransformation2D(matrix=mat_sky2det)


    transform = models.Shift(-output_center[0]) & models.Shift(-output_center[1]) | \
              aff | models.Shift(input_center[0]) & models.Shift(input_center[1])
    return transform


def linear_from_pcf_det2sky(pcffile):
    with open(pcffile) as f:
        lines = [l.strip() for l in f.readlines()]
    factors = lines[lines.index('*Factor 2') + 1].split()
    rotation_angle = float(lines[lines.index('*Rotation') + 1])
    input_rot_center = lines[lines.index('*InputRotationCentre 2') + 1].split()
    output_rot_center = lines[lines.index('*OutputRotationCentre 2') + 1].split()

    det2sky = homothetic_det2sky(input_rot_center, rotation_angle, factors, output_rot_center)
    return det2sky


def fore2asdf(pcffore, outname, ref_kw):
    """
    forward direction : msa 2 ote
    backward_direction: msa 2 fpa

    """
    with open(pcffore) as f:
        lines = [l.strip() for l in f.readlines()]

    fore_det2sky = linear_from_pcf_det2sky(pcffore)
    fore_linear = fore_det2sky
    fore_linear.inverse = fore_det2sky.inverse & Identity(1)

    # compute the polynomial
    degree = int(lines[lines.index('*FitOrder') + 1])

    xcoeff_index = lines.index('*xForwardCoefficients 21 2')
    xlines = lines[xcoeff_index + 1: xcoeff_index + 22]
    xcoeff_forward = coeffs_from_pcf(degree, xlines)
    # Polynomial Correction in x
    x_poly_backward = models.Polynomial2D(degree, name="x_poly_backward", **xcoeff_forward)
    xlines_distortion = lines[xcoeff_index + 22: xcoeff_index + 43]
    xcoeff_forward_distortion = coeffs_from_pcf(degree, xlines_distortion)
    x_poly_backward_distortion = models.Polynomial2D(degree, name="x_backward_distortion",
                                                     **xcoeff_forward_distortion)
    # do chromatic correction
    # the input is Xote, Yote, lam
    model_x_backward = (Mapping((0, 1), n_inputs=3) | x_poly_backward) + \
                     ((Mapping((0,1), n_inputs=3) | x_poly_backward_distortion) * \
                     (Mapping((2,)) | Identity(1)))

    ycoeff_index = lines.index('*yForwardCoefficients 21 2')
    ycoeff_forward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
    y_poly_backward = models.Polynomial2D(degree, name="y_poly_backward",  **ycoeff_forward)

    ylines_distortion = lines[ycoeff_index + 22: ycoeff_index + 43]
    ycoeff_forward_distortion = coeffs_from_pcf(degree, ylines_distortion)
    y_poly_backward_distortion = models.Polynomial2D(degree, name="y_backward_distortion",
                                                     **ycoeff_forward_distortion)

    # do chromatic correction
    # the input is Xote, Yote, lam
    model_y_backward = (Mapping((0,1), n_inputs=3) | y_poly_backward) + \
                     ((Mapping((0, 1), n_inputs=3) | y_poly_backward_distortion) * \
                     (Mapping((2,)) | Identity(1) ))

    xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
    xcoeff_backward = coeffs_from_pcf(degree, lines[xcoeff_index + 1: xcoeff_index + 22])
    x_poly_forward = models.Polynomial2D(degree,name="x_poly_forward", **xcoeff_backward)

    xcoeff_backward_distortion = coeffs_from_pcf(degree, lines[xcoeff_index + 22: xcoeff_index + 43])
    x_poly_forward_distortion = models.Polynomial2D(degree, name="x_forward_distortion", **xcoeff_backward_distortion)

    # the chromatic correction is done here
    # the input is Xmsa, Ymsa, lam
    model_x_forward = (Mapping((0,1), n_inputs=3) | x_poly_forward) + \
                    ((Mapping((0,1), n_inputs=3) | x_poly_forward_distortion) * \
                    (Mapping((2,)) | Identity(1)))

    ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
    ycoeff_backward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
    y_poly_forward = models.Polynomial2D(degree, name="y_poly_forward",**ycoeff_backward)

    ycoeff_backward_distortion = coeffs_from_pcf(degree, lines[ycoeff_index + 22: ycoeff_index + 43])
    y_poly_forward_distortion = models.Polynomial2D(degree, name="y_forward_distortion",
                                                    **ycoeff_backward_distortion)

    # do chromatic correction
    # the input is Xmsa, Ymsa, lam
    model_y_forward = (Mapping((0,1), n_inputs=3) | y_poly_forward) + \
                    ((Mapping((0,1), n_inputs=3) | y_poly_forward_distortion) * \
                    (Mapping((2,)) | Identity(1)))

    #assign inverse transforms
    model_x = model_x_forward.copy()
    model_y = model_y_forward.copy()

    model_x.inverse = model_x_backward
    model_y.inverse = model_y_backward

    output2poly_mapping = Identity(2, name="output_mapping")
    output2poly_mapping.inverse = Mapping([0, 1, 2, 0, 1, 2])
    input2poly_mapping = Mapping([0, 1, 2, 0, 1, 2], name="input_mapping")
    input2poly_mapping.inverse = Identity(2)

    model_poly = input2poly_mapping  | (model_x & model_y) | output2poly_mapping


    model = model_poly | fore_linear

    f = AsdfFile()
    f.tree = ref_kw.copy()
    f.tree['model'] = model
    asdffile = f.write_to(outname)
    return asdffile


def prism2asdf(prifile, tiltyfile, tiltxfile, outname):
    """Create a NIRSPEC prism disperser reference file in ASDF format.

    Combine information stored in disperser_G?.dis and disperser_G?_TiltY.gtp
    files delievred by the IDT.


    Parameters
    ----------
    prifile : list or str
        File with primary information for the PRSIM
    tiltyfile : str
        File with tilt_Y data, e.g. disperser_PRISM_TiltY.gtp.
    tiltxfile: str
        File with tilt_x data, e.g. disperser_PRISM_TiltX.gtp.
    outname : str
        Name of output ASDF file.

    Returns
    -------
    fasdf : asdf.AsdfFile

    """

    params = common_reference_file_keywords("PRISM", "NIRSPEC PRISM Model")
    flist = [prifile, tiltyfile, tiltxfile]

    # translate the files
    for fname in flist:
        try:
            refparams = dict_from_file(fname)
        except:
            print("Disperser file was not converted.")
            raise

        pdict = {}
        coeffs = {}
        parts = fname.lower().split(".")[0]
        ref = str("_".join(parts.split("_")[1:]))

        if "pri" not in fname:
            try:
                for i, c in enumerate(refparams['CoeffsTemperature00']):
                    coeffs['c' + str(i)] = c
                pdict['tilt_model'] = models.Polynomial1D(len(coeffs)-1, **coeffs)
                del refparams['CoeffsTemperature00']
            except KeyError:
                print("Missing CoeffsTemperature in {0}".format(fname))
                raise

        # store the rest of the keys
        for k, v in refparams.items():
            pdict[k] = v
        print(pdict)
        params[ref] = pdict

    fasdf = AsdfFile()
    fasdf.tree = params
    fasdf.write_to(outname)
    return fasdf

def disperser2asdf(disfile, tiltyfile, tiltxfile, outname, ref_kw):
    """
    Create a NIRSPEC disperser reference file in ASDF format.

    Combine information stored in disperser_G?.dis and disperser_G?_TiltY.gtp
    files delievred by the IDT.

    disperser2asdf("disperser_G140H.dis", "disperser_G140H_TiltY.gtp", "disperserG140H.asdf")

    Parameters
    ----------
    disfile : list or str
        A list of .dis files or a wild card string (*.dis).
    tiltyfile : str
        File with tilt_Y data, e.g. disperser_G395H_TiltY.gtp.
    outname : str
        Name of output ASDF file.

    Returns
    -------
    fasdf : asdf.AsdfFile

    """
    disperser = disfile.split('.dis')[0].split('_')[1]
    with open(disfile) as f:
        lines=[l.strip() for l in f.readlines()]
    d = dict.fromkeys(['groove_density', 'theta_z', 'theta_y', 'theta_x', 'tilt_y'])
    d.update(ref_kw)
    for line in lines:
        ind = lines.index('*GRATINGNAME')
        grating_name = lines[ind + 1]
        ind = lines.index('*GROOVEDENSITY')
        d['groove_density'] = float(lines[ind + 1])
        ind = lines.index('*THETAZ')
        d['theta_z'] = float(lines[ind + 1]) / 3600. # in degrees
        ind = lines.index('*THETAX')
        d['theta_x'] = float(lines[ind + 1]) / 3600. # in degrees
        ind = lines.index('*THETAY')
        d['theta_y'] = float(lines[ind + 1]) / 3600. # in degrees
        ind = lines.index('*TILTY')
        d['tilt_y'] = float(lines[ind + 1]) # in degrees
        try:
            ind = lines.index('*TILTX')
            d['tilt_x'] = float(lines[ind + 1]) # in degrees
        except ValueError:
            d['tilt_x'] = 0.0


    assert grating_name in tiltyfile
    assert grating_name in tiltxfile

    with open(tiltyfile) as f:
        s = f.read()
    tiltyd = {}
    vals = s.split('*')
    for line in vals[::-1]:
        if line.startswith("CoeffsTemperature00"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = {}
            for i , c in enumerate([float(c) for c in l[1:1+n]]):
                coeffs['c' + str(i)] = c
            tiltyd['tilt_model'] = models.Polynomial1D(n-1, **coeffs)
        elif line.startswith("Temperatures"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = l[1:1+n]
            tiltyd['temperatures'] = [float(c) for c in coeffs]
        elif line.startswith("Zeroreadings"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = l[1:1+n]
            tiltyd['zeroreadings'] = [float(c) for c in coeffs]
        elif line.startswith("Unit"):
            tiltyd['unit'] = line.split('\n')[1]

    with open(tiltxfile) as f:
        s = f.read()
    tiltxd = {}
    vals = s.split('*')
    for line in vals[::-1]:
        if line.startswith("CoeffsTemperature00"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = {}
            for i , c in enumerate([float(c) for c in l[1:1+n]]):
                coeffs['c' + str(i)] = c
            tiltxd['tilt_model'] = models.Polynomial1D(n-1, **coeffs)
        elif line.startswith("Temperatures"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = l[1:1+n]
            tiltxd['temperatures'] = [float(c) for c in coeffs]
        elif line.startswith("Zeroreadings"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = l[1:1+n]
            tiltxd['zeroreadings'] = [float(c) for c in coeffs]
        elif line.startswith("Unit"):
            tiltxd['unit'] = line.split('\n')[1]


    d['gwa_tiltx'] = tiltyd
    d['gwa_tilty'] = tiltxd
    fasdf = AsdfFile()
    fasdf.tree = d
    fasdf.write_to(outname)
    return fasdf


def wavelength_range(spectral_conf, outname, ref_kw):
    """
    Parameters
    ----------
    spectral_conf : str
        reference file: spectralconfigurations.txt
    outname : str
        output file name
    """
    with open(spectral_conf) as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines][13 :]
    lines = [l.split() for l in lines]
    tree = ref_kw.copy()
    filter_grating = {}
    for l in lines:
        f_g = l[0] + '_' + l[1]
        filter_grating[f_g] = {'order': int(l[2]), 'range': [float(l[3]), float(l[4])]}
    tree['filter_grating'] = filter_grating
    fasdf = AsdfFile()

    fasdf.tree = tree
    fasdf.write_to(outname)
    return fasdf


def msa2asdf(msafile, outname, ref_kw):
    """
    Create an asdf reference file with the MSA description.

    mas2asfdf("MSA.msa", "msa.asdf")

    Parameters
    ----------
    msafile : str
        A fits file with MSA description (MSA.msa)
    outname : str
        Name of output ASDF file.
    """
    f = fits.open(msafile)
    tree = ref_kw.copy()
    data = f[5].data # SLITS and IFU
    header = f[5].header
    shiftx = models.Shift(header['SLITXREF'], name='slit_xref')
    shifty = models.Shift(header['SLITYREF'], name='slit_yref')
    slitrot = models.Rotation2D(header['SLITROT'], name='slit_rot')

    tree[5] = {}
    tree[5]['model'] = slitrot | shiftx & shifty
    tree[5]['data'] = f[5].data
    for i in range(1, 5):
        header = f[i].header
        shiftx = models.Shift(header['QUADXREF'], name='msa_xref')
        shifty = models.Shift(header['QUADYREF'], name='msa_yref')
        slitrot = models.Rotation2D(header['QUADROT'], name='msa_rot')
        tree[i] = {}
        tree[i]['model'] = slitrot | shiftx & shifty
        tree[i]['data'] = f[i].data

    f.close()
    fasdf = AsdfFile()
    fasdf.tree = tree
    fasdf.write_to(outname)
    return fasdf


def fpa2asdf(fpafile, outname, ref_kw):
    """
    Create an asdf reference file with the FPA description.

    The CDP2 delivery includes a fits file - "FPA.fpa" which is the
    input to this function. This file is converted to asdf and is a
    reference file of type "FPA".

    nirspec_fs_ref_tools.fpa2asdf('Ref_Files/CoordTransform/Description/FPA.fpa', 'fpa.asdf')

    Parameters
    ----------
    fpafile : str
        A fits file with FPA description (FPA.fpa)
    outname : str
        Name of output ASDF file.
    """
    with open(fpafile) as f:
        lines = [l.strip() for l in f.readlines()]

    # NRS1
    ind = lines.index("*SCA491_PitchX")
    nrs1_pitchx = float(lines[ind+1])
    ind = lines.index("*SCA491_PitchY")
    nrs1_pitchy = float(lines[ind+1])
    ind = lines.index("*SCA491_RotAngle")
    nrs1_angle = np.rad2deg(float(lines[ind+1]))
    ind = lines.index("*SCA491_PosX")
    nrs1_posx = float(lines[ind+1])
    ind = lines.index("*SCA491_PosY")
    nrs1_posy = float(lines[ind+1])

    # NRS2
    ind = lines.index("*SCA492_PitchX")
    nrs2_pitchx = float(lines[ind+1])
    ind = lines.index("*SCA492_PitchY")
    nrs2_pitchy = float(lines[ind+1])
    ind = lines.index("*SCA492_RotAngle")
    nrs2_angle = np.rad2deg(float(lines[ind+1]))
    ind = lines.index("*SCA492_PosX")
    nrs2_posx = float(lines[ind+1])
    ind = lines.index("*SCA492_PosY")
    nrs2_posy = float(lines[ind+1])

    tree = ref_kw.copy()
    nrs1_sky2det = models.Shift(-nrs1_posx) & models.Shift(-nrs1_posy) | \
                 models.Rotation2D(-nrs1_angle) | \
                 models.Scale(1/nrs1_pitchx) & models.Scale(1/nrs1_pitchy)
    nrs1_det2sky = models.Rotation2D(nrs1_angle) | \
                 models.Scale(nrs1_pitchx) & models.Scale(nrs1_pitchy) | \
                 models.Shift(nrs1_posx) & models.Shift(nrs1_posy)
    nrs1_det2sky.inverse = nrs1_sky2det

    nrs2_sky2det = models.Shift(-nrs2_posx) & models.Shift(-nrs2_posy) | \
                 models.Rotation2D(-nrs2_angle) | \
                 models.Scale(1/nrs2_pitchx) & models.Scale(1/nrs2_pitchy)
    nrs2_det2sky = models.Rotation2D(nrs2_angle) | \
                 models.Scale(nrs2_pitchx) & models.Scale(nrs2_pitchy) | \
                 models.Shift(nrs2_posx) & models.Shift(nrs2_posy)
    nrs2_det2sky.inverse = nrs1_sky2det

    tree['NRS1'] = nrs1_det2sky
    tree['NRS2'] = nrs2_sky2det
    fasdf = AsdfFile()
    fasdf.tree = tree
    fasdf.write_to(outname)
    return fasdf


def ote2asdf(otepcf, outname, ref_kw):
    """
    ref_kw = common_reference_file_keywords('OTE', 'NIRSPEC OTE transform - CDP4')

    ote2asdf('Model/Ref_Files/CoordTransform/OTE.pcf', 'jwst_nirspec_ote_0001.asdf', ref_kw)
    """
    with open(otepcf) as f:
        lines = [l.strip() for l in f.readlines()]

    factors = lines[lines.index('*Factor 2 1') + 1].split()
    # this corresponds to modeling Rotation direction as is
    rotation_angle = float(lines[lines.index('*Rotation') + 1])
    input_rot_center = lines[lines.index('*InputRotationCentre 2 1') + 1].split()
    output_rot_center = lines[lines.index('*OutputRotationCentre 2 1') + 1].split()

    mlinear = homothetic_det2sky(input_rot_center, rotation_angle, factors, output_rot_center)

    degree = int(lines[lines.index('*FitOrder') + 1])

    xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
    xlines = lines[xcoeff_index + 1].split('\t')
    xcoeff_backward = coeffs_from_pcf(degree, xlines)
    x_poly_forward = models.Polynomial2D(degree, name='x_poly_forward', **xcoeff_backward)

    xcoeff_index = lines.index('*xForwardCoefficients 21 2')
    xlines = lines[xcoeff_index + 1].split('\t')
    xcoeff_forward = coeffs_from_pcf(degree, xlines)
    x_poly_backward = models.Polynomial2D(degree, name='x_poly_backward', **xcoeff_forward)

    ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
    ylines = lines[ycoeff_index + 1].split('\t')
    ycoeff_backward = coeffs_from_pcf(degree, ylines)
    y_poly_forward = models.Polynomial2D(degree, name='y_poly_forward', **ycoeff_backward)

    ycoeff_index = lines.index('*yForwardCoefficients 21 2')
    ylines = lines[ycoeff_index + 1].split('\t')
    ycoeff_forward = coeffs_from_pcf(degree, ylines)
    y_poly_backward = models.Polynomial2D(degree, name='y_poly_backward', **ycoeff_forward)

    x_poly_forward.inverse = x_poly_backward
    y_poly_forward.inverse = y_poly_backward

    output2poly_mapping = Identity(2, name='output_mapping')
    output2poly_mapping.inverse = Mapping([0, 1, 0, 1])
    input2poly_mapping = Mapping([0, 1, 0, 1], name='input_mapping')
    input2poly_mapping.inverse = Identity(2)

    model_poly = input2poly_mapping | (x_poly_forward & y_poly_forward) | output2poly_mapping

    model = model_poly | mlinear


    f = AsdfFile()
    f.tree = ref_kw.copy()
    f.tree['model'] = model
    f.write_to(outname)
    return model_poly, mlinear


def ifu_slicer2asdf(ifuslicer, outname):
    """
    Create an asdf reference file with the MSA description.

    ifu_slicer2asdf("IFU_slicer.sgd", "ifu_slicer.asdf")

    Parameters
    ----------
    ifuslicer : str
        A fits file with the IFU slicer description
    outname : str
        Name of output ASDF file.
    """
    ref_kw = common_reference_file_keywords("IFUSLICER", "NIRSPEC IFU SLICER description - CDP4")
    f = fits.open(ifuslicer)
    tree = ref_kw.copy()
    data = f[1].data
    header = f[1].header
    shiftx = models.Shift(header['XREF'], name='ifu_slicer_xref')
    shifty = models.Shift(header['YREF'], name='ifu_slicer_yref')
    rot = models.Rotation2D(header['ROT'], name='ifu_slicer_rot')
    model = rot | shiftx & shifty
    tree['model'] = model
    tree['data'] = f[1].data
    f.close()
    fasdf = AsdfFile()
    fasdf.tree = tree
    fasdf.write_to(outname)
    return fasdf


def ifupost2asdf(ifupost_files, outname):
    """
    Create a reference file of type ``ifupost`` .

    Combines all IDT ``IFU-POST`` reference files in one ASDF file.

    forward direction : MSA to Collimator
    backward_direction: Collimator to MSA

    Parameters
    ----------
    ifupost_files : list
        Names of all ``IFU-POST`` IDT reference files
    outname : str
        Name of output ``ASDF`` file
    """
    ref_kw = common_reference_file_keywords("IFUPOST", "NIRSPEC IFU-POST transforms - CDP4")
    fa = AsdfFile()
    fa.tree = ref_kw
    for fifu in ifupost_files:
        n = int((fifu.split('IFU-POST_')[1]).split('.pcf')[0])
        fa.tree[n] = {}
        with open(fifu) as f:
            lines = [l.strip() for l in f.readlines()]
        factors = lines[lines.index('*Factor 2') + 1].split()
        rotation_angle = float(lines[lines.index('*Rotation') + 1])
        input_rot_center = lines[lines.index('*InputRotationCentre 2') + 1].split()
        output_rot_center = lines[lines.index('*OutputRotationCentre 2') + 1].split()
        linear_sky2det = homothetic_sky2det(input_rot_center, rotation_angle, factors, output_rot_center)

        degree = int(lines[lines.index('*FitOrder') + 1])

        xcoeff_index = lines.index('*xForwardCoefficients 21 2')
        xlines = lines[xcoeff_index + 1: xcoeff_index + 22]
        xcoeff_forward = coeffs_from_pcf(degree, xlines)
        x_poly_forward = models.Polynomial2D(degree, name='x_poly_forward', **xcoeff_forward)

        ycoeff_index = lines.index('*yForwardCoefficients 21 2')
        ycoeff_forward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
        y_poly_forward = models.Polynomial2D(degree, name='y_poly_forward', **ycoeff_forward)

        xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
        xcoeff_backward = coeffs_from_pcf(degree, lines[xcoeff_index + 1: xcoeff_index + 22])
        x_poly_backward = models.Polynomial2D(degree, name='x_poly_backward', **xcoeff_backward)

        ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
        ycoeff_backward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
        y_poly_backward = models.Polynomial2D(degree, name='y_poly_backward', **ycoeff_backward)

        output2poly_mapping = Identity(2, name='output_mapping')
        output2poly_mapping.inverse = Mapping([0, 1, 0, 1])
        input2poly_mapping = Mapping([0, 1, 0, 1], name='input_mapping')
        input2poly_mapping.inverse = Identity(2)

        model_poly = input2poly_mapping | (x_poly_forward & y_poly_forward) | output2poly_mapping

        model = linear_sky2det | model_poly
        fa.tree[n]['model'] = model
    asdffile = fa.write_to(outname)
    return asdffile


ref_files = "/internal/1/astropy/embray-compound-mdboom-gwcs3/models/nirspec-model-2014"


def nirspec_models_to_asdf():
    ref = {}
    camera_refname = os.path.join(ref_files, "CoordTransform", "Camera.pcf")
    camera_name = "jwst_nirspec_camera_0001.asdf"
    ref_kw = common_reference_file_keywords("CAMERA", "NIRSPEC Camera Model - CDP4")
    try:
        pcf2asdf(camera_refname, camera_name, ref_kw)
    except:
        print("Camera file was not converted.")
        raise
    ref["CAMERA"] = camera_name

    ref_kw = common_reference_file_keywords("COLLIMATOR", "NIRSPEC Collimator Model - CDP4")
    collimator_refname = os.path.join(ref_files, "CoordTransform", "Collimator.pcf")
    collimator_name = "jwst_nirspec_collimator_0001.asdf"
    try:
        pcf2asdf(collimator_refname, collimator_name, ref_kw)
    except:
        print("Collimator file was not converted.")
        raise
    ref["COLLIMATOR"] = collimator_name

    ref_kw = common_reference_file_keywords("WAVELENGTHRANGE", "NIRSPEC Spectral Configurations - CDP4")
    wavelengthrange_name = "jwst_nirspec_wavelengthrange_0001.asdf"
    wavelength_range('spectralconfigurations.txt', wavelengthrange_name, ref_kw)
    ref["WAVELENGTHRANGE"] = wavelengthrange_name

    ref_kw = common_reference_file_keywords("MSA", "NIRSPEC MSA Description - CDP4")
    mas_refname = os.path.join(ref_files, "Description", "MSA.msa")
    msa_name = "jwst_nirspec_msa_0001.asdf"
    try:
        msa2asdf(mas_refname, msa_name, ref_kw)
    except:
        print("MSA file was not converted")
        raise
    ref['MSA'] = msa_name

    ref_kw = common_reference_file_keywords("FPA", "NIRSPEC FPA Description - CDP4")
    fpa_refname = os.path.join(ref_files, "Description", "FPA.fpa")
    fpa_name = "jwst_nirspec_fpa_0001.asdf"
    try:
        fpa2asdf(fpa_refname, fpa_name, ref_kw)
    except:
        print("FPA file was not converted.")
        raise
    ref["FPA"] = fpa_name

    for i, filter in enumerate(["CLEAR", "F070LP", "F100LP", "F110W", "F140X", "F170LP", "F290LP"]):
        ref_kw = common_reference_file_keywords("FORE", "NIRSPEC FORE Model - CDP4")
        ref_kw['filter'] = filter
        filename = "Fore_{0}.pcf".format(filter)
        fore_refname = os.path.join(ref_files, "CoordTransform", filename)
        fore_name = "jwst_nirspec_fore_000{0}.asdf".format(str(i+1))
        try:
            fore2asdf(fore_refname, fore_name, ref_kw)
        except:
            print(("FORE file was not created - filter {0}".format(filter)))
            raise
        ref[filter] = fore_name

    for i, grating in enumerate(["G140H", "G140M", "G235H", "G235M", "G395H", "G395M", "MIRROR"]):
        # disperser refers to the whole reference file, either on disk or asdf
        # grating to the parameter the reference file is selected on
        ref_kw = common_reference_file_keywords("DISPERSER", "NIRSPEC Disperser Description - CDP4")
        ref_kw['grating'] = grating
        dis_file = "disperser_" + grating + ".dis"
        gtpy_file = "disperser_" + grating + "_TiltY.gtp"
        gtpx_file = "disperser_" + grating + "_TiltX.gtp"
        disp_refname = os.path.join(ref_files, "Description", dis_file)
        tilty_refname = os.path.join(ref_files, "Description", gtpy_file)
        tiltx_refname = os.path.join(ref_files, "Description", gtpx_file)
        disperser_name =  "jwst_nirspec_disperser_000{0}.asdf".format(str(i+1))

        try:
            disperser2asdf(disp_refname, tilty_refname, tiltx_refname, disperser_name, ref_kw)
        except:
            print("Disperser file was not converted.")
            raise
        ref["DISPERSER"] = disperser_name

    ref_kw = common_reference_file_keywords("OTE", "NIRSPEC OTE Model - CDP4")
    ote_refname = os.path.join(ref_files, "CoordTransform", "OTE.pcf")
    ote_name = "jwst_nirspec_ote_0001.asdf"
    try:
        ote2asdf(ote_refname, ote_name, ref_kw)
    except:
        print("OTE file was not converted.")
        raise
    ref["OTE"] = ote_name
