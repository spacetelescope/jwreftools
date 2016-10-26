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


"""
from asdf import AsdfFile
from astropy.modeling.models import Mapping

import read_siaf_table

def create_nircam_distortion(detector, aperture, outname):
    """
    Create an asdf reference file with all distortion components for the NIRCam imager.

    NOTE: The IDT has not provided any distortion information. The files are constructed
    using ISIM transformations provided/(computed?) by the TEL team which they use to
    create the SIAF file.
    These reference files should be replaced when/if the IDT provides us with distortion.

    Parameters
    ----------
    detector : str
        NRCB1, NRCB2, NRCB3, NRCB4, NRCB5, NRCA1, NRCA2, NRCA3, NRCA4, NRCA5
        (NRCA5 and NRCB5 are the LW detectors) 
    aperture : str
        Name of the aperture/subarray. (e.g. FULL, SUB160, SUB320, SUB640, GRISM_F322W2)
    outname : str
        Name of output file.

    Examples
    --------

    """
    numdet = detector[-1]
    module = detector[-2]
    channel = 'SW'
    if numdet == '5':
        channel = 'LW'

    full_aperture = detector + '_' + aperture


    #"Forward' transformations. science --> ideal --> V2V3
    sci2idlx, sci2idly, sciunit, idlunit = read_siaf_table.get_siaf_transform(full_aperture,'science','ideal', 5)
    idl2v2v3x, idl2v2v3y = read_siaf_table.get_siaf_v2v3_transform(full_aperture,from_system='ideal')

    #'Reverse' transformations. V2V3 --> ideal --> science
    v2v32idlx, v2v32idly = read_siaf_table.get_siaf_v2v3_transform(full_aperture,to_system='ideal')
    idl2scix, idl2sciy, idlunit, sciunit = read_siaf_table.get_siaf_transform(full_aperture,'ideal','science', 5)
    
 
    #Map the models together to make a single transformation
    model =  Mapping([0, 1, 0, 1]) | sci2idlx & sci2idly | Mapping([0, 1, 0, 1]) | idl2v2v3x & idl2v2v3y
    model_inv =  Mapping([0, 1, 0, 1]) | v2v32idlx & v2v32idly | Mapping([0, 1, 0, 1]) | idl2scix & idl2sciy
    model.inverse = model_inv

    tree = {"title": "NIRCAM Distortion",
            "instrument": "NIRCAM",
            "pedigree": "GROUND",
            "reftype" : "DISTORTION",
            "author": "B. Hilbert",
            "detector": detector,
            "module": module,
            "channel": channel,
            "subarray": aperture,
            "description": "Distortion model function created from SIAF coefficients",
            "exp_type": "NRC_IMAGE",
            "useafter": "2014-01-01T00:00:00",
            "model": model
            }

    fasdf = AsdfFile()
    fasdf.tree = tree
    fasdf.write_to(outname)
