"""
This module contains functions to create NIRCAM reference files.

NIRCAM model:
im.meta.instrument.name :'NIRCAM'
im.meta.instrument.channel : 'SHORT'
im.meta.instrument.module : 'B'
im.meta.instrument.detector : 'NRCB1'
im.meta.exposure.type : 'NRC_IMAGE'

???
im.meta.instrument.pupil : 'FLAT' (would be GRISMR or GRISMC for slitless)

Transform Paths for Imaging mode:

science --> ideal --> V2V3
V2V3 --> ideal --> science

Where the "science" frame has units of distorted pixels. The "ideal" frame
is distortion-free and is the distance in arcseconds from 0,0. V2/V3 is also
in units of arcseconds from the refrence pixel.

Grism
-----
Each NIRCam module has 2 grisms that disperse spectra in orthogonal directions,
along detector rows (GrismR) and columns (GrismC).

For each observation, a grism in the pupil wheel is used in combination with
a filter in the  filter wheel.   Broader filters yield longer spectra
(in detector pixels). The available filters depend on the observing mode


"""
import re
import numpy as np
import datetime
from asdf import AsdfFile
from asdf.tags.core import Software, HistoryEntry

from astropy.io import fits
from astropy import units as u
from astropy.modeling.models import Polynomial1D

import read_siaf_table
from jwst.datamodels import NIRCAMGrismModel
from jwst.datamodels import wcs_ref_models


def common_reference_file_keywords(reftype, title, description, exp_type,
                                   useafter, author, module=None, fname=None,
                                   pupil=None, **kwargs):
    """
    exp_type can be also "N/A", or "ANY".
    """
    ref_file_common_keywords = {
        "AUTHOR": author,
        "DESCRIP": description,
        "EXP_TYPE": exp_type,
        "INSTRUME": "NIRCAM",
        "PEDIGREE": "GROUND",
        "REFTYPE": reftype,
        "TELESCOP": "JWST",
        "TITLE": title,
        "USEAFTER": useafter,
        }


    if fname is not None:
        ref_file_common_keywords["FILTER"] = fname
    if pupil is not None:
        ref_file_common_keywords["PUPIL"] = pupil
    if module is not None:
        ref_file_common_keywords["MODULE"] = module

    ref_file_common_keywords.update(kwargs)
    return ref_file_common_keywords


def create_grism_config(conffile="",
                        pupil=None,
                        module=None,
                        author="STScI",
                        history="",
                        outname=""):
    """
    Create an asdf reference file to hold Grism C (column) or Grism R (rows)
    configuration, no sensativity information is included

    Note: The orders are named alphabetically, i.e. Order A, Order B
    There are also sensativity fits files which are tables of wavelength,
    sensativity, and error. These are specified in the conffile but will
    not be read in and saved in the output reference file for now.
    It's possible they may be included in the future, either here or as
    a separate reference files. Their use here would be to help define the
    min and max wavelengths which set the extent of the dispersed trace on
    the grism image. Convolving the sensitiviy file with the filter throughput
    allows one to calculate the wavelength of minimum throughput which defines
    the edges of the trace.

    direct_filter is not specified because it assumes that the wedge
    information (wx,wy) is included in the conf file in one of the key-value
    pairs, where the key includes the beam designation

     this reference file also contains the polynomial model which is appropriate
     for the coefficients which are listed.
     wavelength = DISPL(order,x0,y0,t)
     dx = DISPX(order,x0,y0,t)
     dy = DISPY(order,x0,y0,t)

     t = INVDISPX(order,x0,y0,dx)
     t = INVDISPY(order,x0,y0,dy)
     t = INVDISL(order,x0,y0, wavelength)



    Parameters
    ----------
    conffile : str
        The text file with configuration information
    pupil : str
        Name of the grism the conffile corresponds to
    module : str
        Name of the module
    author : str
        The name of the author
    history : str
        A comment about the refrence file to be saved with the meta information
    outname : str
        Output name for the reference file

    Returns
    -------
    fasdf : asdf.AsdfFile(jwst.datamodels.NIRCAMGrismModel)

    """
    if not history:
        history = "Created from {0:s}".format(conffile)

    # if pupil is none get from filename like NIRCAM_modB_R.conf
    if pupil is None:
        pupil = "GRISM" + conffile.split(".")[0][-1]
    # if module is none get from filename
    if module is None:
        module = conffile.split(".")[0][-3]
    print("Pupil is {}".format(pupil))

    ref_kw = common_reference_file_keywords("specwcs",
                                            "NIRCAM Grism Parameters",
                                            "{0:s} dispersion models".format(pupil),
                                            "NRC_GRISM",
                                            "2014-01-01T00:00:00",
                                            author,
                                            model_type="NIRCAMGrismModel",
                                            module=module,
                                            pupil=pupil)

    # get all the key-value pairs from the input file
    conf = dict_from_file(conffile)
    beamdict = split_order_info(conf)
    # beam = re.compile('^(?:[+\-]){0,1}[a-zA-Z0-9]{0,1}$')  # match beam only
    # read in the sensitivity tables to save their content
    # they currently have names like this: NIRCam.A.1st.sensitivity.fits
    # translated as inst.beam/order.param
    temp = dict()
    etoken = re.compile("^[a-zA-Z]*_(?:[+\-]){1,1}[1,2]{1,1}")  # find beam key
    for b, bdict in beamdict.items():
            temp[b] = dict()
            # for key in bdict:
                # if 'SENSITIVITY' in key:
                #     print("Reading sensitivity from: {0:s}".format(bdict[key]))
                #     sdata = read_sensitivity_file(bdict[key])
                #     skeys = sdata.keys()
                #     for k in skeys:
                #         temp[b][key+"_"+k] = sdata[k]

    # add the new beam information to beamdict and remove spurious beam info
    for k in temp:
        for kk in temp[k]:
            if etoken.match(kk):
                kk = kk.replace("_{}".format(k), "")
            beamdict[k][kk] = temp[k][kk]

    # add min and max mag info if not provided
    # also make beam coeff lists
    # wx are the wedge offsets for the filters

    # for k, bdict in beamdict.items():
    #     if isinstance(bdict, dict):
    #         keys = bdict.keys()
    #         if "MMAG_EXTRACT" not in keys:
    #             beamdict[k]["MMAG_EXTRACT"] = 99.0
            # if maxmag not in keys:
            #    beamdict[k][maxmag] = 0.0
            # if "wx" not in keys:
            #    beamdict[k]['wx'] = 0.0
            # if "wy" not in keys:
            #    beamdict[k]['wy'] = 0.0
            # add the model for transforms



    # for NIRCAM, the R and C grism coefficients contain zeros where
    # the dispersion is in the opposite direction. Meaning, the GRISMR,
    # which disperses along ROWS has coefficients of zero in the y models
    # and vice versa.
    #
    # There are separate reference files for each grism. Depending on the grism
    # dispersion direction you either want to use the dx from source center or
    # the dy from source center in the inverse dispersion relationship which is
    # used to calculate the t value needed to calculate the wavelength at that pixel.
    #
    # The model creation here takes all of this into account by looking at the GRISM[R/C]
    # the file is used for and creating a reference model with the appropriate dispersion
    # direction in use. This eliminates having to decide which direction to calculate
    # the dispersion from given the input x,y pixel in the dispersed image.
    orders = beamdict.keys()

    # dispersion models valid per order and direction saved to reference file
    # Forward
    invdispl = []
    invdispx = []
    invdispy = []
    # Backward
    displ = []
    dispx = []
    dispy = []

    for order in orders:
        # convert the displ wavelengths to microns if the input file is still in angstroms
        l0 = beamdict[order]['DISPL'][0] / 10000.
        l1 = beamdict[order]['DISPL'][1] / 10000.

        # create polynomials using the coefficients of each order

        # This holds the wavelength lookup coeffs
        # This model is  INVDISPL for backward and returns t
        # This model should be DISPL for forward and returns wavelength
        if l1 == 0:
            lmodel = Polynomial1D(1, c0=0, c1=0)
        else:
            lmodel = Polynomial1D(1, c0=-l0/l1, c1=1./l1)
        invdispl.append(lmodel)
        lmodel = Polynomial1D(1, c0=l0, c1=l1)
        displ.append(lmodel)

        # This holds the x coefficients, for the R grism this model is the
        # the INVDISPX returning t, for the C grism this model is the DISPX
        x0, x1 = beamdict[order]['DISPX']
        xmodel = Polynomial1D(1, c0=x0, c1=x1)
        dispx.append(xmodel)
        if x1 == 0:
            xmodel = Polynomial1D(1, c0=0, c1=0)
        else:
            xmodel = Polynomial1D(1, c0=-x0/x1, c1=1./x1)
        invdispx.append(xmodel)

        # This holds the y coefficients, for the C grism, this model is
        # the INVDISPY, returning t, for the R grism, this model is the DISPY
        y0, y1 = beamdict[order]['DISPY']
        ymodel = Polynomial1D(1, c0=y0, c1=y1)
        dispy.append(ymodel)
        if y1 == 0:
            ymodel = Polynomial1D(1, c0=0, c1=0)
        else:
            ymodel = Polynomial1D(1, c0=-y0/y1, c1=1./y1)
        invdispy.append(ymodel)

    # change the orders into translatable integers
    # so that we can look up the order with the proper index
    oo = [int(o) for o in beamdict]

    ref = NIRCAMGrismModel()
    ref.meta.update(ref_kw)
    ref.meta.input_units = u.micron
    ref.meta.output_units = u.micron
    ref.displ = displ
    ref.dispx = dispx
    ref.dispy = dispy
    ref.invdispx = invdispx
    ref.invdispy = invdispy
    ref.invdispl = invdispl
    ref.orders = oo
    entry = HistoryEntry({'description': history, 'time': datetime.datetime.utcnow()})
    sdict = Software({'name': 'nircam_reftools.py',
             'author': author,
             'homepage': 'https://github.com/spacetelescope/jwreftools',
             'version': '0.7.1'})
    entry['sofware'] = sdict
    ref.history = [entry]
    ref.to_asdf(outname)
    ref.validate()


def create_grism_waverange(outname="",
                           history="Ground NIRCAM Grismwavelengthrange",
                           author="STScI",
                           filter_range=None):
    """Create a wavelengthrange reference file.

    Supply a filter range dictionary keyed on order or use the default

    """
    ref_kw = common_reference_file_keywords("wavelengthrange",
                                            "NIRCAM Grism wavelenghtrange",
                                            "NIRCAM Grism+Filter Wavelength Ranges",
                                            "NRC_GRISM",
                                            "2014-01-01T00:00:00",
                                            author,
                                            model_type="WavelengthrangeModel",
                                            )

    if filter_range is None:
        # These numbers from Nor Pirzkal, in microns
        filter_range = {1: {'F250M': [2.500411072, 4.800260833],
                            'F277W': [2.500411072, 3.807062006],
                            'F300M': [2.684896869, 4.025318456],
                            'F322W2': [2.5011293930000003, 4.215842089],
                            'F335M': [3.01459734, 4.260432726],
                            'F356W': [3.001085025, 4.302320901],
                            'F360M': [3.178096344, 4.00099629],
                            'F410M': [3.6267051809999997, 4.5644598],
                            'F430M': [4.04828939, 4.511761774],
                            'F444W': [3.696969216, 4.899565197],
                            'F460M': [3.103778615, 4.881999188],
                            'F480M': [4.5158154679999996, 4.899565197]},
                        2: {'F250M': [2.500411072, 2.667345336],
                            'F277W': [2.500411072, 3.2642254050000004],
                            'F300M': [2.6659796289999997, 3.2997071729999994],
                            'F322W2': [2.5011293930000003, 4.136119434],
                            'F335M': [2.54572003, 3.6780519760000003],
                            'F356W': [2.529505253, 4.133416971],
                            'F360M': [2.557881113, 4.83740855],
                            'F410M': [2.5186954019999996, 4.759037127],
                            'F430M': [2.5362614100000003, 4.541488865],
                            'F444W': [2.5011293930000003, 4.899565197],
                            'F460M': [2.575447122, 4.883350419],
                            'F480M': [2.549773725, 4.899565197]}}

    # array of integers
    orders = list(filter_range.keys())
    orders.sort()

    # same filters for every order, array of strings
    wrange_selector = list(filter_range[orders[0]].keys())
    wrange_selector.sort()

    # The lists below need
    # to remain ordered to be correctly referenced
    wavelengthrange = []
    for order in orders:
        o = []
        for fname in wrange_selector:
            o.append(filter_range[order][fname])
        wavelengthrange.append(o)

    ref = wcs_ref_models.WavelengthrangeModel()
    ref.meta.update(ref_kw)
    ref.meta.input_units = u.micron
    ref.meta.output_units = u.micron
    ref.wrange_selector = wrange_selector
    ref.wrange = wavelengthrange
    ref.order = orders

    entry = HistoryEntry({'description': history, 'time': datetime.datetime.utcnow()})
    sdict = Software({'name': 'nircam_reftools.py',
                      'author': author,
                      'homepage': 'https://github.com/spacetelescope/jwreftools',
                      'version': '0.7.1'})
    entry['sofware'] = sdict
    ref.history = [entry]
    ref.to_asdf(outname)
    ref.validate()



def create_filter_transmission(filename="",
                             outname="",
                             filtername="",
                             detector="",
                             history="",
                             author=""):
    """Read in a filter transmission file and save a passband reference file.
    Assumes that wavelength is in microns

    filename: str
        Name of the text file with wavelength + transmission columns

    outname: str
        The name of the output asdf reference file

    filtername: str
        The filter this data represents

    history: str
        History comment for the reference file

    author: str
        Name of the author


    Notes:
    ------
    The passband file is used with the sensativity information for modeling
    grism observations. It's not currently used for extraction or with the
    dispersion polynomials
    """

    wave, transmission = np.loadtxt(filename, skiprows=1,
                                    usecols=(0, 1), unpack=True)
    # wave *= 10000.  # convert to angstroms

    # Put these in the reference file?
    numdet = detector[-1]
    module = detector[-2]
    channel = 'SHORT'
    if numdet == '5':
        channel = 'LONG'
        detector = detector[0:4] + 'LONG'

    # In the reference file headers, we need to switch NRCA5 to NRCALONG,
    # and the same for module B.

    if not history:
        history = "Created from {0:s}".format(filename)

    tree = {"TITLE": "NIRCAM Passband Transmission",
            "TELESCOP": "JWST",
            "INSTRUME": "NIRCAM",
            "PEDIGREE": "GROUND",
            "REFTYPE": "TRANSMISSION",
            "AUTHOR": author,
            "DETECTOR": detector,
            "MODULE": module,
            "CHANNEL": channel,
            "SUBARRAY": opsgname,
            "DESCRIP": "Filter Transmission for {0:s}".format(filtername),
            "EXP_TYPE": "NRC_IMAGE",
            "FILTER": filtername,
            "USEAFTER": "2014-01-01T00:00:00",
            "WAVELENGTH_UNITS": u.micron,
            }
    tree["wavelength"] = wave
    tree["transmission"] = transmission

    fasdf = AsdfFile()
    fasdf.tree = tree
    sdict = {'name': 'nircam_reftools.py',
             'author': author,
             'homepage': 'https://github.com/spacetelescope/jwreftools',
             'version': '0.7'}

    fasdf.add_history_entry(history, software=sdict)
    fasdf.write_to(outname)


def read_sensitivity_file(filename):
    """
    Read the sensativity fits file.
    This is assumed to be in the same format Nor gave me which is an MEF
    file with a table in the first extension. Also assumes that there is
    a column called WAVELENGTH and SENSITIVITY.

    Parameters
    ----------
    filename : str
        The name of the fits file

    Returns
    -------
    A dictionary of lists where the keys are the column names

    Notes
    -----
    The sensitivity files contain the grism sensitivie measures by wavelength.
    The can be convolved with the filter throughput files in order to find
    the minimum wavelength and maximum wavelength for which throughput flux
    falls to zero. These wavelengths are then used to define the edge extents
    of the dispersed trace on the grism image for each object.
    """
    if not isinstance(filename, str):
        raise ValueError("Expected the name of the sensitivity fits file")

    sens = dict()
    with fits.open(filename) as fh:
        columns = fh[1].header['TTYPE*']
        for c in columns.values():
            sens[c] = fh[1].data.field(c)

    # pin the edges to zero
    sens['SENSITIVITY'][0:2] = 0.
    sens['SENSITIVITY'][-2:] = 0.

    sens["WRANGE"] = [np.min(sens['WAVELENGTH'][sens['SENSITIVITY'] != 0]),
                      np.max(sens['WAVELENGTH'][sens['SENSITIVITY'] != 0])]
    return sens


def create_nircam_distortion(coefffile, detector, aperture, opgsname, outname):
    """Create an asdf reference file with all distortion components for the NIRCam imager.

    NOTE: The IDT has not provided any distortion information. The files are
    constructed using ISIM transformations provided/(computed?) by the TEL team
    which they use to create the SIAF file. These reference files should be
    replaced when/if the IDT provides us with distortion.

    Parameters
    ----------
    detector : str
        NRCB1, NRCB2, NRCB3, NRCB4, NRCB5, NRCA1, NRCA2, NRCA3, NRCA4, NRCA5
    aperture : str
        Name of the aperture/subarray. (e.g. FULL, SUB160, SUB320,
        SUB640, GRISM_F322W2)
    outname : str
        Name of output file.

    Examples
    --------

    """
    numdet = detector[-1]
    module = detector[-2]
    channel = 'SHORT'
    if numdet == '5':
        channel = 'LONG'

    full_aperture = detector + '_' + aperture

    # "Forward' transformations. science --> ideal --> V2V3
    sci2idlx, sci2idly, sciunit, idlunit = read_siaf_table.get_siaf_transform(coefffile,full_aperture,'science','ideal', 5)
    idl2v2v3x, idl2v2v3y = read_siaf_table.get_siaf_v2v3_transform(coefffile,full_aperture,from_system='ideal')

    # 'Reverse' transformations. V2V3 --> ideal --> science
    v2v32idlx, v2v32idly = read_siaf_table.get_siaf_v2v3_transform(coefffile,full_aperture,to_system='ideal')
    idl2scix, idl2sciy, idlunit, sciunit = read_siaf_table.get_siaf_transform(coefffile,full_aperture,'ideal','science', 5)


    # Map the models together to make a single transformation
    model =  Mapping([0, 1, 0, 1]) | sci2idlx & sci2idly | Mapping([0, 1, 0, 1]) | idl2v2v3x & idl2v2v3y
    model_inv =  Mapping([0, 1, 0, 1]) | v2v32idlx & v2v32idly | Mapping([0, 1, 0, 1]) | idl2scix & idl2sciy
    model.inverse = model_inv


    # In the reference file headers, we need to switch NRCA5 to NRCALONG, and same
    # for module B.
    if detector[-1] == '5':
        detector = detector[0:4] + 'LONG'

    meta = {"TITLE": "NIRCAM Distortion",
            "TELESCOP": "JWST",
            "INSTRUME": "NIRCAM",
            "PEDIGREE": "GROUND",
            "REFTYPE": "DISTORTION",
            "AUTHOR": "B. Hilbert",
            "DETECTOR": detector,
            "MODULE": module,
            "CHANNEL": channel,
            "SUBARRAY": opgsname,
            "DESCRIP": "Distortion model function created from SIAF coefficients",
            "EXP_TYPE": "NRC_IMAGE",
            "USEAFTER": "2014-01-01T00:00:00",
            "model": model
            }

    fasdf = AsdfFile()
    fasdf.tree = tree

    sdict = {'name': 'nircam_reftools.py',
             'author': 'B.Hilbert',
             'homepage': 'https://github.com/spacetelescope/jwreftools',
             'version': '0.7'}

    fasdf.add_history_entry("File created from a file of distortion "
                            "coefficients, NIRCam_SIAF_2016-09-29.csv, "
                            "provided by Colin Cox in October 2016. "
                            "Software used: https://github.com/spacetele"
                            "scope/jwreftools", software=sdict)

    fasdf.write_to(outname)


def split_order_info(keydict):
    """
    Designed to take as input the dictionary created by dict_from_file and for
    nircam, split out and accumulate the keys for each beam/order.
    The keys must have the beam in their string, the spurious beam designation
    is removed from the returned dictionary. Keywords with the same first name
    in the underscore separated string followed by a number are assumed to be
    ranges


    Parameters
    ----------
    keydict : dictionary
        Dictionary of key value pairs

    Returns
    -------
    dictionary of beams, where each beam has a dictionary of key-value pairs
    Any key pairs which are not associated with a beam get a separate entry
    """

    if not isinstance(keydict, dict):
        raise ValueError("Expected an input dictionary")

    # has beam name fits token
    token = re.compile('^[a-zA-Z]*_(?:[+\-]){0,1}[a-zA-Z0-9]{0,1}_*')
    rangekey = re.compile('^[a-zA-Z]*_[0-1]{1,1}$')
    rdict = dict()  # return dictionary
    beams = list()

    # prefetch number of Beams, beam is the second string
    for key in keydict:
        if token.match(key):
            b = key.split("_")[1].upper()
            if b not in beams:
                beams.append(b)
    for b in beams:
        rdict[b] = dict()

    #  assumes that keys are sep with underscore and beam is in second section
    for key in keydict:
        if not token.match(key):
            rdict[key] = keydict[key]  # not associated with a beam
        else:
            b = key.split("_")[1].upper()
            newkey = key.replace("_{}".format(b), "")
            rdict[b][newkey] = keydict[key]

    # look for range variables to make them into tuples
    for b, d in rdict.items():
        keys = d.keys()
        rkeys = []
        odict = {}
        for k in keys:
            if rangekey.match(k):
                rkeys.append(k)
        for k in rkeys:
            mlist = [m for m in rkeys if k.split("_")[0] in m]
            root = mlist[0].split("_")[0]
            if root not in odict:
                for mk in mlist:
                    if eval(mk[-1]) == 0:
                        zero = d[mk]
                    elif eval(mk[-1]) == 1:
                        one = d[mk]
                    else:
                        raise ValueError("Unexpected range variable {}"
                                         .format(mk))
                odict[root] = (zero, one)
        # combine the dictionaries and remove the old keys
        d.update(odict)
        for k in rkeys:
            del d[k]

    return rdict


def dict_from_file(filename):
    """Read in a file and return a named tuple of the key value pairs.

    This is a generic read for a text file with the line following format:

    keyword<token>value

    Where keyword should start with a character, not a number
    Non-alphabetic starting characters are ignored
    <token> can be space or comma

    Parameters
    ----------
    filename : str
        Name of the file to interpret

    Examples
    --------
    dict_from_file('NIRCAM_C.conf')

    Returns
    -------
    dictionary of deciphered keys and values

    """
    token = '\s+|(?<!\d)[,](?!\d)'
    letters = re.compile("(^[a-zA-Z])")  # starts with a letter
    numbers = re.compile("(^(?:[+\-])?(?:\d*)(?:\.)?(?:\d*)?(?:[eE][+\-]?\d*$)?)")
    empty = re.compile("(^\s*$)")  # is a blank line

    print("\nReading {0:s}  ...".format(filename))
    with open(filename, 'r') as fh:
        lines = fh.readlines()
    content = dict()
    for line in lines:
        value = None
        key = None
        if not empty.match(line):
            if letters.match(line):
                pair = re.split(token, line.strip(), maxsplit=3)
                if len(pair) == 2:  # key and value exist
                    key = pair[0]  # first item is the key
                    val = pair[1]  # second item is the value
                    if letters.match(val):
                        value = val
                    if numbers.fullmatch(val):
                        value = eval(val)
                if len(pair) == 3:  # key min max exist
                    key = pair[0]
                    val1, val2 = pair[1:]
                    if numbers.fullmatch(val1) and numbers.fullmatch(val2):
                        value = (eval(val1), eval(val2))
                    else:
                        raise ValueError("Min/max values expected for {0}"
                                         .format(key))
        # ignore the filter file pointings and the sensitivity files these are
        # used for simulation
        if key and (value is not None):
            if (("FILTER" not in key) and ("SENSITIVITY" not in key)):
                content[key] = value
                print("Setting {0:s} = {1}".format(key, value))

    return content
