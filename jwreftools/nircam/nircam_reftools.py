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

