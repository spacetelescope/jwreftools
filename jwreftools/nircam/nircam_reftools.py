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
from astropy.modeling.models import Mapping,Shift
from astropy.io import ascii
import read_siaf_table
import numpy as np

def create_nircam_distortion(coefffile, detector, aperture, exp_type, pupil, outname):
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
    aperture : str
        Name of the aperture/subarray. (e.g. FULL, SUB160, SUB320, SUB640, GRISM_F322W2)
    outname : str
        Name of output file.

    Examples
    --------

    """
    #Read in the SIAF table
    siaf = ascii.read(coefffile,header_start=1)
    
    numdet = detector[-1]
    module = detector[-2]
    channel = 'SHORT'
    if numdet == '5':
        channel = 'LONG'
        
    full_aperture = detector + '_' + aperture

    #Find the row that matches the requested aperture
    match = siaf['AperName'] == full_aperture
    if np.any(match) == False:
        print("Aperture name {} not found in input CSV file.".format(aperture))
        sys.exit()

    siaf_row = siaf[match]
    
    #"Forward' transformations. science --> ideal --> V2V3

    #Frist, find the offset between 0,0 and the reference location, in pixels
    xshift,yshift = read_siaf_table.get_refpix(siaf_row)

    #Now the core of the model. 2D-Polynomial to go from units of science pixels
    #to the 'ideal' (distortion-free) coordinate system
    sci2idlx, sci2idly, sciunit, idlunit = read_siaf_table.get_siaf_transform(siaf_row,full_aperture,'science','ideal', 5)

    #And another 2D polynomial to go from 'ideal' coordinates to V2,V3
    idl2v2v3x, idl2v2v3y = read_siaf_table.get_siaf_v2v3_transform(siaf_row,full_aperture,from_system='ideal')

    #Finally, we need to shift by the v2,v3 value of the reference location in order to get to
    #absolute v2,v3 coordinates
    v2shift,v3shift = read_siaf_table.get_v2v3ref(siaf_row)
    
    #'Reverse' transformations. V2V3 --> ideal --> science
    #In this case we just need the inverses of the polynomials from above
    #Inverses of the shifts are known by astropy modeling 
    v2v32idlx, v2v32idly = read_siaf_table.get_siaf_v2v3_transform(siaf_row,full_aperture,to_system='ideal')
    idl2scix, idl2sciy, idlunit, sciunit = read_siaf_table.get_siaf_transform(siaf_row,full_aperture,'ideal','science', 5)


    #Now create a compound model for each with the appropriate inverse
    sci2idl = Mapping([0,1,0,1]) | sci2idlx & sci2idly
    sci2idl.inverse = Mapping([0,1,0,1]) | idl2scix & idl2sciy

    idl2v2v3 = Mapping([0,1,0,1]) | idl2v2v3x & idl2v2v3y
    idl2v2v3.inverse = Mapping([0,1,0,1]) | v2v32idlx & v2v32idly
    
    #Now string the models together to make a single transformation
    
    #First the core of the translations, without the shifts
    #Scale the models by 1/3600 because the SIAF has things in units of arcsec
    #but the pipeline assumes degrees

    #Official word is that reffiles should now work in arcsec,
    #not degrees, so remove the Scale(1./3600), but we also need to
    #account for the difference of 1 between the SIAF coordinate values
    #(indexed to 1) and python (indexed to 0). Nadia said that this
    #shift should be present in the distortion reference file.

    core_model = sci2idl | idl2v2v3 #| Scale(1./3600) & Scale(1./3600)
    
    #Now add in the shifts to create the full model
    #including the shift to go from 0-indexed python coords to 1-indexed
    #SIAF coords
    index_shift = Shift(1)
    model = index_shift & index_shift | xshift & yshift | core_model | v2shift & v3shift

    #Since the inverse of all model components are now defined, the total model
    #inverse is also defined automatically

    #core_model_inv =  Mapping([0, 1, 0, 1]) | v2v32idlx & v2v32idly | Mapping([0, 1, 0, 1]) | idl2scix & idl2sciy | Scale(1./3600) & Scale(1./3600)
    #Now add in the inverse shifts
    #model_inv = v2shift.inverse & v3shift.inverse | core_model_inv | xshift.inverse & yshift.inverse
    #model.inverse = model_inv


    #In the reference file headers, we need to switch NRCA5 to NRCALONG, and same
    #for module B.
    if detector[-1] == '5':
        detector = detector[0:4] + 'LONG'

    tree = {"TITLE": "NIRCAM Distortion",
            "TELESCOP": "JWST",
            "INSTRUMENT": "NIRCAM",
            "PEDIGREE": "GROUND",
            "REFTYPE" : "DISTORTION",
            "AUTHOR": "B. Hilbert",
            "LITREF": "https://github.com/spacetelescope/jwreftools",
            "DETECTOR": detector,
            "MODULE": module,
            "CHANNEL": channel,
            "PUPIL": pupil,
            "SUBARRAY": 'ANY',
            "DESCRIP": "Distortion model function created from SIAF coefficients",
            "EXP_TYPE": exp_type,
            "USEAFTER": "2014-01-01T00:00:00",
            "model": model
            }

    fasdf = AsdfFile()
    fasdf.tree = tree

    sdict = {'name':'nircam_reftools.py','author':'B.Hilbert','homepage':'https://github.com/spacetelescope/jwreftools','version':'1.0'}
    
    fasdf.add_history_entry("File created from distortion coefficients contained in SIAF ("+coefffile+"), provided by Colin Cox in October 2016. Software used: https://github.com/spacetelescope/jwreftools",software=sdict)

    fasdf.write_to(outname)
