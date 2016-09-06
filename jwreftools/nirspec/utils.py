import numpy as np
from astropy.modeling import models
from astropy.modeling.models import Mapping, Identity
from asdf import AsdfFile


def homothetic_det2sky(input_center, angle, scale, output_center, name=""):
    """
    Create the homothetic transform from a .pcf file.

    The forward direction is sky to detector.
    Parameters
    ----------
    input_center : ndarray of shape (1, 2) or iterable
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

    rotmat_det2sky = [[np.cos(angle), -np.sin(angle)],
                      [np.sin(angle), np.cos(angle)]]
    scaling = np.array([[1/scale[0] , 0], [0, 1/scale[1]]])
    mat_det2sky = np.dot(rotmat_det2sky, scaling)

    aff = models.AffineTransformation2D(matrix=mat_det2sky, name="affine_det2sky_{0}".format(name))


    transform = models.Shift(-output_center[0], name="x_input_center_det2sky") & models.Shift(-output_center[1],  name="y_input_center_det2sky") | \
              aff | models.Shift(input_center[0],  name="x_output_center_det2sky") & models.Shift(input_center[1],  name="y_output_center_det2sky")
    return transform


def homothetic_sky2det(input_center, angle, scale, output_center, name=""):
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

    aff = models.AffineTransformation2D(matrix=mat_sky2det, name="affine_{0}".format(name))


    transform = models.Shift(-input_center[0], name="x_input_center_{0}".format(name)) & models.Shift(-input_center[1], name="y_input_center_{0}".format(name)) | aff | \
              models.Shift(output_center[0], name="x_output_center_{0}".format(name)) & models.Shift(output_center[1], name="y_output_center_{0}".format(name))
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
    f.add_history_entry("Build 6")
    f.tree['model'] = model
    f.write_to(outname)


def common_reference_file_keywords(reftype, title, description, exp_type,
                                   useafter, author, filename, **kwargs):
    """
    exp_type can be also "N/A", or "ANY".
    """
    ref_file_common_keywords = {
        "author": author,
        "description": description,
        "exp_type": exp_type,
        "filename": filename,
        "instrume": "NIRSPEC",
        "pedigree": "GROUND",
        "reftype": reftype,
        "telescope": "JWST",
        "title": title,
        "useafter": useafter,
        }
    ref_file_common_keywords.update(kwargs)
    return ref_file_common_keywords
