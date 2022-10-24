
import re
import datetime

import numpy as np

from asdf.tags.core import Software, HistoryEntry

from astropy import units as u
from astropy.modeling.models import Polynomial1D, Polynomial2D, Const2D

from jwst.datamodels import NIRCAMGrismModel
from jwst.datamodels import wcs_ref_models


def common_reference_file_keywords(reftype=None,
                                   title=None,
                                   description=None,
                                   exp_type=None,
                                   author="STScI",
                                   useafter="2014-01-01T00:00:00",
                                   module=None,
                                   fname=None,
                                   pupil=None, **kwargs):
    """
    exp_type can be also "N/A", or "ANY".
    """
    if exp_type is None:
        raise ValueError("exp_type not set")
    if reftype is None:
        raise ValueError("Expected reftype value")

    ref_file_common_keywords = {
        "author": author,
        "description": description,
        "exposure": {"type": exp_type},
        "instrument": {"name": "NIRCAM"},
        "pedigree": "ground",
        "reftype": reftype,
        "telescope": "JWST",
        "title": title,
        "useafter": useafter,
        }

    if fname is not None:
        ref_file_common_keywords["instrument"]["filter"] = fname
    if pupil is not None:
        ref_file_common_keywords["instrument"]["pupil"] = pupil
    if module is not None:
        ref_file_common_keywords["instrument"]["module"] = module

    ref_file_common_keywords.update(kwargs)
    return ref_file_common_keywords


def dict_from_file(conffile):
    """Read in the CONF file into a dictionary."""

    def update_beamdict(beamdict, beam, dispx, dispy, displ):
        print(f"updating beam {beam}")
        beamdict[beam] = {'DISPX': np.array([line.split()[1:] for line in dispx], dtype=float),
                          'DISPY': np.array([line.split()[1:] for line in dispy], dtype=float),
                          'DISPL': np.array([line.split()[1:] for line in displ], dtype=float)
                          }
    f = open(conffile)
    lines = f.readlines()
    f.close()
    beamdict = {}
    beam = ""
    for line in lines:
        if line.startswith('BEAM'):

            if beam.strip():
                update_beamdict(beamdict, beam, dispx, dispy, displ)
            dispx = []
            dispy = []
            displ = []
            beam = line.split('_')[1].strip()

        elif line.startswith('DISPL'):
            displ.append(line)
        elif line.startswith('DISPX'):
            dispx.append(line)
        elif line.startswith('DISPY'):
            dispy.append(line)
        else:
            continue

    # Capture the last order reading from the CONF file
    update_beamdict(beamdict, beam, dispx, dispy, displ)
    return beamdict


def get_degree_coeffs(e_grismconf):
    """ Return the Polynomial degree and a dict of coefficients.

    Parameters
    ----------
    e-grismconf : list
        A list of numbers - one line with coefficients from the CONF file.
    """

    # keys are the size of the list of coefficients, values - the polynomial degree.
    size_degree_map = {1: 0,
                       6: 2,
                       10: 3
                       }
    degree = size_degree_map[e_grismconf.size]
    e0 = e_grismconf
    if e_grismconf.size == 1:
        coeffs = {'c0_0': e_grismconf[0]}
    elif e_grismconf.size == 10:
        coeffs = {'c0_0': e_grismconf[0], 'c1_0': e_grismconf[1], 'c2_0': e_grismconf[3], 'c3_0': e_grismconf[6],
                  'c0_1': e_grismconf[2], 'c1_1': e_grismconf[4], 'c0_2': e_grismconf[5], 'c0_3': e_grismconf[9],
                  'c2_1': e_grismconf[7], 'c1_2': e_grismconf[8]}
    return size_degree_map[e_grismconf.size], coeffs


def get_polynomial(e_grismconf):
    """Return an astropy.modeling Polynomial from coefficients in the CONF file.

    Parameters
    ----------
    e-grismconf : list
        A list of numbers - one line with coefficients from the CONF file.

    Returns
    -------
    poly : `~astropy.modeling.polynomial.Polynomial`
        A Polynomial object.
    """
    degree, coeffs = get_degree_coeffs(e_grismconf)
    return Polynomial2D(degree, **coeffs)


def create_grism_config(conffile,
                        filter_name,
                        pupil,
                        author="STScI",
                        history="",
                        outname="test.asdf"):
    """

    Parameters
    ----------
    conffile : str
        The NIRCam CONF file, aXe format
    filter_name : str
        The name of the filter
    pupil : str
        The grism name, one of [GRISMC, GRISMR]
    author : str
        The name of the author.
    history : str
        A comment about the refrence file to be saved with
        the meta information.
    outname : str
        The name of the output ASDF file.

    Create a reference file of reftype "specwcs" to store the transforms
    from direct image to grism image. The file is in ASDF format.
    No sensitivity information is included.
    The original CONF files are in https://github.com/npirzkal/GRISM_NIRCAM.
    There are also sensitivity files in FITS format. These are converted to
    "photom" reference files by the NIRCam team.


    Returns
    -------
    fasdf : asdf.AsdfFile(jwst.datamodels.NIRCamGrismModel)

    Examples
    --------
    >>> create_grism_config('NIRCAM_F360M_modA_C.conf', 'F360M', 'GRISMC', author="STScI",
                             outname="test_nircam_specwcs.asdf")
    """

    if not history:
        history = "Created from {0:s}".format(conffile)

    ref_kw = common_reference_file_keywords(reftype="specwcs",
                                            description="{0:s} dispersion model parameters".format(pupil),
                                            exp_type="NRC_WFSS",
                                            model_type='NIRCAMGrismModel',
                                            pupil=pupil,
                                            filtername=filter_name,
                                            history=history,
                                            author=author,
                                            filename=outname,
                                            )

    # get all the key-value pairs from the input file
    beamdict = dict_from_file(conffile)

    # The lists below need
    # to remain ordered and referenced by filter or order
    orders = sorted(beamdict.keys())
    displ = []
    dispx = []
    dispy = []
    invdispl = []
    int_orders = [int(order) for order in orders]

    """
    V4 of the NIRCam reference files have updated transforms for beam "+1"
    from commisioning data. Beam "+2" has not been calibrated in flight yet
    and stores transforms from gound testing. The transforms for beam "+1"
    (order 1) have field dependence of the traces as well as the wavelength
    solutions.
    """
    for order in orders:
        e = beamdict[order]['DISPL']
        if len(e) == 3:
            e0, e1, e2 = e
            cpoly_0 = get_polynomial(e0)
            cpoly_1 = get_polynomial(e1)
            cpoly_2 = get_polynomial(e2)
            displ.append((cpoly_0, cpoly_1, cpoly_2))
        elif len(e) == 2:
            e0, e1 = e
            cpoly_0 = Polynomial1D(1, c0=e0[0], c1=e1[0])
            displ.append((cpoly_0,))

        e = beamdict[order]['DISPX']

        if len(e) == 3:
            # final poly is ans = cpoly_0 + t*cpoly_1 + t**2 * cpoly_2
            e0, e1, e2 = e
            cpoly_0 = get_polynomial(e0)
            cpoly_1 = get_polynomial(e1)
            cpoly_2 = get_polynomial(e2)

            dispx.append((cpoly_0, cpoly_1, cpoly_2))
        elif len(e) == 2:
            e0, e1 = e
            cpoly_0 = Polynomial1D(1, c0=e0[0], c1=e1[0])
            dispx.append((cpoly_0,))

        e = beamdict[order]['DISPY']
        if len(e) == 3:
            e0, e1, e2 = e
            cpoly_0 = get_polynomial(e0)
            cpoly_1 = get_polynomial(e1)
            cpoly_2 = get_polynomial(e2)
            dispy.append((cpoly_0, cpoly_1, cpoly_2))
        elif len(e) == 2:
            e0, e1 = e
            cpoly_0 = Polynomial1D(1, c0=e0[0], c1=e1[0])
            dispy.append((cpoly_0,))

    ref = NIRCAMGrismModel()
    ref.meta.instance.update(ref_kw)
    ref.meta.exposure.p_exptype = "NRC_WFSS|NRC_TSGRISM"
    ref.meta.input_units = u.micron
    ref.meta.output_units = u.micron
    ref.meta.input_units = u.micron
    ref.meta.output_units = u.micron
    ref.dispx = dispx
    ref.dispy = dispy
    ref.displ = displ
    ref.invdispl = invdispl

    ref.orders = [int(order) for order in orders]
    entry = HistoryEntry({'description': history,
                          'time': datetime.datetime.utcnow()})
    sdict = Software({'name': 'nircam_grism_reffiles.py',
                      'author': "STScI",
                      'homepage': 'https://github.com/spacetelescope/jwreftools',
                      'version': '0.7.1'})
    entry['software'] = sdict
    ref.history.append(entry)
    ref.to_asdf(outname)
    ref.validate()
    return ref


def create_tsgrism_wavelengthrange(outname="nircam_tsgrism_wavelengthrange.asdf",
                                   history="Ground NIRCAM TSGrism wavelengthrange",
                                   author="STScI",
                                   wavelengthrange=None,
                                   extract_orders=None):
    """Create a wavelengthrange reference file for NIRCAM TSGRISM mode.

    Parameters
    ----------
    outname: str
        The output name of the file
    history: str
        History information about it's creation
    author: str
        Person or entity making the file
    wavelengthrange: list(tuples)
        A list of tuples that set the order, filter, and
        wavelength range min and max
    extract_orders: list[list]
        A list of lists that specify

    """
    ref_kw = common_reference_file_keywords(reftype="wavelengthrange",
                                            title="NIRCAM TSGRISM reference file",
                                            description="NIRCAM Grism-Filter Wavelength Ranges",
                                            exp_type="NRC_TSGRISM",
                                            author=author,
                                            pupil="ANY",
                                            model_type="WavelengthrangeModel",
                                            filename=outname,
                                            )

    if wavelengthrange is None:
        # This is a list of tuples that specify the
        # order, filter, wave min, wave max
        wavelengthrange = [(1, 'F277W', 2.500411072, 3.807062006),
                           (1, 'F322W2', 2.5011293930000003, 4.215842089),
                           (1, 'F356W', 3.001085025, 4.302320901),
                           (1, 'F444W', 3.696969216, 4.899565197),
                           (2, 'F277W', 2.500411072, 3.2642254050000004),
                           (2, 'F322W2', 2.5011293930000003, 4.136119434),
                           (2, 'F356W', 2.529505253, 4.133416971),
                           (2, 'F444W', 2.5011293930000003, 4.899565197),
                           ]

    # array of integers of unique orders
    orders = sorted(set((x[0] for x in wavelengthrange)))
    filters = sorted(set((x[1] for x in wavelengthrange)))

    # Nircam has not specified any limitation on the orders
    # that should be extracted by default yet so all are
    # included.
    if extract_orders is None:
        extract_orders = [('F277W', [1]),
                          ('F322W2', [1]),
                          ('F356W', [1]),
                          ('F444W', [1]),
                          ]

    ref = wcs_ref_models.WavelengthrangeModel()
    ref.meta.update(ref_kw)
    ref.meta.exposure.p_exptype = "NRC_TSGRISM"
    ref.meta.input_units = u.micron
    ref.meta.output_units = u.micron
    ref.wavelengthrange = wavelengthrange
    ref.extract_orders = extract_orders
    ref.order = orders
    ref.waverange_selector = filters

    history = HistoryEntry({'description': history,
                            'time': datetime.datetime.utcnow()})
    software = Software({'name': 'nircam_reftools.py',
                         'author': author,
                         'homepage': 'https://github.com/spacetelescope/jwreftools',
                         'version': '0.7.1'})
    history['software'] = software
    ref.history = [history]
    ref.validate()
    ref.to_asdf(outname)


def create_wfss_wavelengthrange(outname="nircam_wfss_wavelengthrange.asdf",
                                history="Ground NIRCAM Grism wavelengthrange",
                                author="STScI",
                                wavelengthrange=None,
                                extract_orders=None):
    """Create a wavelengthrange reference file for NIRCAM.

    Parameters
    ----------
    outname: str
        The output name of the file
    history: str
        History information about it's creation
    author: str
        Person or entity making the file
    wavelengthrange: list(tuples)
        A list of tuples that set the order, filter, and
        wavelength range min and max
    extract_orders: list[list]
        A list of lists that specify

    """
    ref_kw = common_reference_file_keywords(reftype="wavelengthrange",
                                            title="NIRCAM WFSS reference file",
                                            description="NIRCAM Grism-Filter Wavelength Ranges",
                                            exp_type="NRC_WFSS",
                                            author=author,
                                            pupil="ANY",
                                            model_type="WavelengthrangeModel",
                                            filename=outname,
                                            )

    if wavelengthrange is None:
        # This is a list of tuples that specify the
        # order, filter, wave min, wave max
        wavelengthrange = [(1, 'F250M', 2.500411072, 4.800260833),
                           (1, 'F277W', 2.500411072, 3.807062006),
                           (1, 'F300M', 2.684896869, 4.025318456),
                           (1, 'F322W2', 2.5011293930000003, 4.215842089),
                           (1, 'F335M', 3.01459734, 4.260432726),
                           (1, 'F356W', 3.001085025, 4.302320901),
                           (1, 'F360M', 3.178096344, 4.00099629),
                           (1, 'F410M', 3.6267051809999997, 4.5644598),
                           (1, 'F430M', 4.04828939, 4.511761774),
                           (1, 'F444W', 3.696969216, 4.899565197),
                           (1, 'F460M', 3.103778615, 4.881999188),
                           (1, 'F480M', 4.5158154679999996, 4.899565197),
                           (2, 'F250M', 2.500411072, 2.667345336),
                           (2, 'F277W', 2.500411072, 3.2642254050000004),
                           (2, 'F300M', 2.6659796289999997, 3.2997071729999994),
                           (2, 'F322W2', 2.5011293930000003, 4.136119434),
                           (2, 'F335M', 2.54572003, 3.6780519760000003),
                           (2, 'F356W', 2.529505253, 4.133416971),
                           (2, 'F360M', 2.557881113, 4.83740855),
                           (2, 'F410M', 2.5186954019999996, 4.759037127),
                           (2, 'F430M', 2.5362614100000003, 4.541488865),
                           (2, 'F444W', 2.5011293930000003, 4.899565197),
                           (2, 'F460M', 2.575447122, 4.883350419),
                           (2, 'F480M', 2.549773725, 4.899565197),
                           ]
    # array of integers of unique orders
    orders = sorted(set((x[0] for x in wavelengthrange)))
    filters = sorted(set((x[1] for x in wavelengthrange)))

    # Nircam has not specified any limitation on the orders
    # that should be extracted by default yet so all are
    # included.
    if extract_orders is None:
        extract_orders = []
        for f in filters:
            extract_orders.append([f, orders])

    ref = wcs_ref_models.WavelengthrangeModel()
    ref.meta.update(ref_kw)
    ref.meta.exposure.p_exptype = "NRC_WFSS"
    ref.meta.input_units = u.micron
    ref.meta.output_units = u.micron
    ref.wavelengthrange = wavelengthrange
    ref.extract_orders = extract_orders
    ref.order = orders
    ref.waverange_selector = filters

    history = HistoryEntry({'description': history,
                            'time': datetime.datetime.utcnow()})
    software = Software({'name': 'nircam_reftools.py',
                         'author': author,
                         'homepage': 'https://github.com/spacetelescope/jwreftools',
                         'version': '0.7.1'})
    history['software'] = software
    ref.history = [history]
    ref.validate()
    ref.to_asdf(outname)


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
        if token.match(key):
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

'''
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
'''
