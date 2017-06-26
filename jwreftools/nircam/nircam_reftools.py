#! /usr/bin/env python
"""
This module contains functions to create NIRCAM distortion reference 
files using the DistortionModel model.

NIRCAM model:
im.meta.instrument.name :'NIRCAM'
im.meta.instrument.channel : 'SHORT'
im.meta.instrument.module : 'B'
im.meta.instrument.detector : 'NRCB1'
im.meta.exposure.type : 'NRC_IMAGE'

im.meta.instrument.pupil : 'FLAT' (would be GRISMR or GRISMC for 
slitless)

Transform Paths for Imaging mode:

science --> ideal --> V2V3 
V2V3 --> ideal --> science

Where the "science" frame has units of distorted pixels. The "ideal" 
frame is distortion-free and is the distance in arcseconds from 0,0. 
V2/V3 is also in units of arcseconds from the refrence pixel.


"""
from asdf import AsdfFile
import numpy as np
import sys
from astropy.modeling.models import Mapping,Shift
from astropy.io import ascii
from astropy import units as u
from jwst.datamodels import DistortionModel,util
import read_siaf_table


def create_nircam_distortion(coefffile, detector, aperture,
                             outname, sci_pupil,
                             sci_subarr, sci_exptype):
    """
    Create an asdf reference file with all distortion components 
    for the NIRCam imager.

    NOTE: The IDT has not provided any distortion information. The 
    files are constructed using ISIM transformations provided by 
    the TEL team which they use to create the SIAF file.

    Parameters
    ----------
    detector : str
        NRCB1, NRCB2, NRCB3, NRCB4, NRCB5, NRCA1, NRCA2, NRCA3, 
        NRCA4, NRCA5
    aperture : str
        Name of the aperture/subarray. (e.g. FULL, SUB160, SUB320, 
        SUB640, GRISM_F322W2)
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
        print("Aperture name {} not found in input CSV file.".
              format(aperture))
        sys.exit()

    siaf_row = siaf[match]
    
    #"Forward' transformations. science --> ideal --> V2V3

    #Frist, find the offset between 0,0 and the reference
    #location, in pixels
    xshift,yshift = read_siaf_table.get_refpix(siaf_row)

    #Now the core of the model. 2D-Polynomial to go from units
    #of science pixels to the 'ideal' (distortion-free)
    #coordinate system
    sci2idlx, sci2idly, sciunit, idlunit = read_siaf_table.get_siaf_transform(siaf_row,full_aperture,'science','ideal', 5)

    #And another 2D polynomial to go from 'ideal' coordinates
    #to V2,V3
    idl2v2v3x, idl2v2v3y = read_siaf_table.get_siaf_v2v3_transform(siaf_row,full_aperture,from_system='ideal')

    #Finally, we need to shift by the v2,v3 value of the reference
    #location in order to get to absolute v2,v3 coordinates
    v2shift,v3shift = read_siaf_table.get_v2v3ref(siaf_row)
    
    #'Reverse' transformations. V2V3 --> ideal --> science
    #In this case we just need the inverses of the polynomials
    #from above
    #Inverses of the shifts are known by astropy modeling 
    v2v32idlx, v2v32idly = read_siaf_table.get_siaf_v2v3_transform(siaf_row,full_aperture,to_system='ideal')
    idl2scix, idl2sciy, idlunit, sciunit = read_siaf_table.get_siaf_transform(siaf_row,full_aperture,'ideal','science', 5)

    #Now create a compound model for each with the appropriate
    #inverse
    sci2idl = Mapping([0,1,0,1]) | sci2idlx & sci2idly
    sci2idl.inverse = Mapping([0,1,0,1]) | idl2scix & idl2sciy

    idl2v2v3 = Mapping([0,1,0,1]) | idl2v2v3x & idl2v2v3y
    idl2v2v3.inverse = Mapping([0,1,0,1]) | v2v32idlx & v2v32idly
    
    #Now string the models together to make a single transformation
    
    #First the core of the translations, without the shifts
    #Scale the models by 1/3600 because the SIAF has things in
    #units of arcsec but the pipeline assumes degrees

    #Official word is that reffiles should now work in arcsec,
    #not degrees, so remove the Scale(1./3600), but we also need
    #to account for the difference of 1 between the SIAF
    #coordinate values (indexed to 1) and python (indexed to 0).
    #Nadia said that this shift should be present in the
    #distortion reference file.

    core_model = sci2idl | idl2v2v3 
    
    #Now add in the shifts to create the full model
    #including the shift to go from 0-indexed python coords to
    #1-indexed
    
    #SIAF coords
    index_shift = Shift(1)
    model = index_shift & index_shift | xshift & yshift | core_model | v2shift & v3shift

    #Since the inverse of all model components are now defined,
    #the total model inverse is also defined automatically

    #In the reference file headers, we need to switch NRCA5 to
    #NRCALONG, and same for module B.
    if detector[-1] == '5':
        detector = detector[0:4] + 'LONG'

    #Save using the DistortionModel datamodel
    d = DistortionModel(model=model, input_units=u.pix,
                        output_units=u.arcsec)

    #Populate metadata

    #keyword values in science data to which this file should
    #be applied
    p_pupil = ''
    for p in sci_pupil:
        p_pupil = p_pupil + p + '|'

    p_subarr = ''
    for p in sci_subarr:
        p_subarr = p_subarr + p + '|'
            
    p_exptype = ''
    for p in sci_exptype:
        p_exptype = p_exptype + p + '|'

    d.meta.instrument.p_pupil = p_pupil
    d.meta.subarray.p_subarray = p_subarr
    d.meta.exposure.p_exptype = p_exptype
    
    #d.meta.instrument.p_pupil = "CLEAR|F162M|F164N|F323N|F405N|F470N|"
    #d.meta.p_subarray = "FULL|SUB64P|SUB160|SUB160P|SUB320|SUB400P|SUB640|SUB32TATS|SUB32TATSGRISM|SUB8FP1A|SUB8FP1B|SUB96DHSPILA|SUB96DHSPILB|SUB64FP1A|SUB64FP1B|"
    #d.meta.exposure.p_exptype = "NRC_IMAGE|NRC_TSIMAGE|NRC_FLAT|NRC_LED|NRC_WFSC|"

    #metadata describing the reference file itself
    d.meta.title = "NIRCam Distortion"
    d.meta.instrument.name = "NIRCAM"
    d.meta.instrument.module = module
    d.meta.instrument.channel = channel
    d.meta.instrument.detector = detector
    d.meta.telescope = 'JWST'
    d.meta.subarray.name = 'FULL'
    d.meta.pedigree = 'GROUND'
    d.meta.reftype = 'DISTORTION' 
    d.meta.author = 'B. Hilbert'
    d.meta.litref = "https://github.com/spacetelescope/jwreftools"
    d.meta.description = "Distortion model from SIAF coefficients in {}".format(coefffile)
    #d.meta.exp_type = exp_type
    d.meta.useafter = "2014-10-01T00:00:00"

    #to be ready for the future where we will have filter-dependent solutions
    d.meta.instrument.filter = 'N/A'

    
    #Create initial HISTORY ENTRY
    
    sdict = {'name':'nircam_reftools.py','author':'B.Hilbert','homepage':'https://github.com/spacetelescope/jwreftools','version':'0.71'}
    
    entry = util.create_history_entry("File created from distortion coefficients contained in SIAF ("+coefffile+"), provided by Colin Cox in March 2017. The format of this distortion reference file matches that used by version 7.1 of the JWST Calibration Pipeline. The translation model converts from units of pixels on the detector to V2,V3 in units of degrees, as well as the inverse. P_PUPIL values have been updated in these files to reflect the correct pupil wheel options.",software=sdict)
    d.history = [entry]


    #Create additional HISTORY entries
    #entry2 = util.create_history_entry("The format of this distortion reference file matches that used by version 7.1 of the JWST Calibration Pipeline. The translation model converts from units of pixels on the detector to V2,V3 in units of degrees, as well as the inverse.")
    #d.history.append(entry2)
    

    d.save(outname)
    print("Output saved to {}".format(outname))



def create_nircam_pixel_area_map(coefffile, detector, aperture,
                             outname, sci_pupil,
                             sci_subarr, sci_exptype):
    """
    Create an asdf reference file with all distortion components 
    for the NIRCam imager.

    NOTE: The IDT has not provided any distortion information. The 
    files are constructed using ISIM transformations provided by 
    the TEL team which they use to create the SIAF file.

    Parameters
    ----------
    detector : str
        NRCB1, NRCB2, NRCB3, NRCB4, NRCB5, NRCA1, NRCA2, NRCA3, 
        NRCA4, NRCA5
    aperture : str
        Name of the aperture/subarray. (e.g. FULL, SUB160, SUB320, 
        SUB640, GRISM_F322W2)
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
        print("Aperture name {} not found in input CSV file.".
              format(aperture))
        sys.exit()

    siaf_row = siaf[match]
    
    #"Forward' transformations. science --> ideal --> V2V3

    #Frist, find the offset between 0,0 and the reference
    #location, in pixels
    xshift,yshift = read_siaf_table.get_refpix(siaf_row)

    #Now the core of the model. 2D-Polynomial to go from units
    #of science pixels to the 'ideal' (distortion-free)
    #coordinate system
    sci2idlx, sci2idly, sciunit, idlunit = read_siaf_table.get_siaf_transform(siaf_row,full_aperture,'science','ideal', 5)

    #And another 2D polynomial to go from 'ideal' coordinates
    #to V2,V3
    idl2v2v3x, idl2v2v3y = read_siaf_table.get_siaf_v2v3_transform(siaf_row,full_aperture,from_system='ideal')

    #Finally, we need to shift by the v2,v3 value of the reference
    #location in order to get to absolute v2,v3 coordinates
    v2shift,v3shift = read_siaf_table.get_v2v3ref(siaf_row)
    
    #'Reverse' transformations. V2V3 --> ideal --> science
    #In this case we just need the inverses of the polynomials
    #from above
    #Inverses of the shifts are known by astropy modeling 
    v2v32idlx, v2v32idly = read_siaf_table.get_siaf_v2v3_transform(siaf_row,full_aperture,to_system='ideal')
    idl2scix, idl2sciy, idlunit, sciunit = read_siaf_table.get_siaf_transform(siaf_row,full_aperture,'ideal','science', 5)

    #Now create a compound model for each with the appropriate
    #inverse
    sci2idl = Mapping([0,1,0,1]) | sci2idlx & sci2idly
    sci2idl.inverse = Mapping([0,1,0,1]) | idl2scix & idl2sciy
