"""This script creates distortion reference files for NIRCam using the coefficients from
pysiaf.
"""
import nircam_distortion_reffiles_from_pysiaf as ref
import numpy as np
import os

nrc_a_apertures = ['NRCA{}_FULL'.format(i+1) for i in range(5)]
nrc_b_apertures = ['NRCB{}_FULL'.format(i+1) for i in range(5)]
nrc_apertures = nrc_a_apertures + nrc_b_apertures

hist = ("This reference file was created from the distortion coefficients contained in pysiaf "
        "(and therefore the PRD) as of 24 Oct 2019. This includes the fix for the round trip (V2,V3 -> x,y -> V2, V3) error "
        "present in previous versions of the coefficients. This update is described in JIRA issues: "
        "JWSTSIAF-161, JWSTSIAF-123. pysiaf version 0.6.1 was used to access the appropriate "
        "distortion coefficients. The translation model converts from units of pixels on the detector to "
        "V2,V3 in units of arcseconds, as well as the inverse.")

basedir = 'reference_files/distortion/24Oct2019_round_trip_error_fixed'

# IMAGING metadata-----------------------------------------------
sw_imaging_pupil = ['CLEAR', 'F162M', 'F164N', 'GDHS0', 'GDHS60', 'WLM8', 'WLP8', 'PINHOLES', 'MASKIPR', 'FLAT']
lw_imaging_pupil = ['CLEAR', 'F323N', 'F405N', 'F466N', 'F470N', 'PINHOLES', 'MASKIPR', 'GRISMR', 'GRISMC', 'FLAT']

# Put grism/wfss exp_types in SW list of exposure types because simultaneous SW data were
# tagged as e.g. NRC_GRISM when grism observations were made with the LW detector
# early on.
sw_exptype = ['NRC_IMAGE', 'NRC_TSIMAGE', 'NRC_FLAT', 'NRC_LED',
              'NRC_WFSC', 'NRC_TACQ', 'NRC_TACONFIRM', 'NRC_FOCUS',
              'NRC_DARK', 'NRC_WFSS', 'NRC_TSGRISM', 'NRC_GRISM']
lw_exptype = sw_exptype #+ ['NRC_WFSS', 'NRC_TSGRISM', 'NRC_GRISM']

for aperture in nrc_apertures:
    detector, apname = aperture.split('_')
    outname = os.path.join(basedir, '{}_distortion.asdf'.format(aperture))

    if np.int(detector[-1]) < 5:
        pupil = sw_imaging_pupil
        #subarr = sw_imaging_subarr
        # Leaving subarrays as any means that new subarrays can be added later and the file won't need to be updated
        subarr = ['GENERIC']
        exp_type = sw_exptype
    else:
        pupil = lw_imaging_pupil
        #subarr = lw_imaging_subarr
        subarr = ['GENERIC']
        exp_type = lw_exptype

    ref.create_nircam_distortion(detector, apname, outname, pupil, subarr, exp_type, hist)

# CORONAGRAPHY metadata------------------------------------------
# For coronagraphy, we expect the bar and round wedges to have different astrometric solutions,
# in the future if not now, so treat these separately. Also, there is currently no aperture in SIAF for
# A1 and A3 that include the wedge, and so to make reference files for them, I would need to manually
# combine the NRCA1_FULL coefficients with the wedge offsets in the SIAF file. John is putting in a
# ticket to add these apertures to SIAF. It might be worth waiting for that to happen before creating
# these reference files?

a2_round_aperture = ['NRCA2_FULL_MASK210R']
a4_bar_aperture = ['NRCA4_FULL_MASKSWB']
a5_round_aperture = ['NRCA5_FULL_MASK335R']
a5_bar_aperture = ['NRCA5_FULL_MASKLWB']

b1_round_aperture = ['NRCB1_MASK210R']
b3_bar_aperture = ['NRCB3_MASKSWB']
b5_round_aperture = ['NRCB5_MASK335R']
b5_bar_aperture = ['NRCB5_MASKLWB']

nrc_coron_apertures = a2_round_aperture + a4_bar_aperture + a5_round_aperture + a5_bar_aperture + b1_round_aperture + b3_bar_aperture + b5_round_aperture + b5_bar_aperture

coron_exptype = ['NRC_CORON', 'NRC_TACQ', 'NRC_TACONFIRM']

for aperture in nrc_coron_apertures:
    parts = aperture.split('_')
    detector = parts[0]
    if len(parts) == 2:
        apname = parts[1]
    elif len(parts) == 3:
        apname = '{}_{}'.format(parts[1], parts[2])

    if aperture[-1] == 'R':
        pupil = ['MASKRND']
    elif aperture[-1] == 'B':
        pupil = ['MASKBAR']
    else:
        raise ValueError('Bad pupil value: {}'.format(aperture))

    outname = os.path.join(basedir, '{}_{}_distortion.asdf'.format(detector, pupil[0]))

    # Leaving subarrays as any means that new subarrays can be added later and the file won't need to be updated
    subarr = ['GENERIC']
    exp_type = coron_exptype

    ref.create_nircam_distortion(detector, apname, outname, pupil, subarr, exp_type, hist)
