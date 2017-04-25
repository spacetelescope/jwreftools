import nircam_reftools as ref
from astropy.io import ascii
import numpy as np

siaf_file = 'NIRCam_SIAF-OSS_Summary_20161117_MMA.csv'
nrc_siaf = 'NIRCam_SIAF_2016-11-21.csv'

siaf = ascii.read(siaf_file,data_start=18)
siaf['col4'].fill_value = 'na'
siaf['col18'].fill_value = 'na'
siaf = siaf.filled()

#Apertures to use for the new case where we only need
#full frame files
apertures = []
mods = ['A','B']
for mod in mods:
    for i in range(5):
        apertures.append('NRC'+mod+str(i+1)+'_FULL')

coron = ['NRCA2_FULL_MASK210R','NRCA5_FULL_MASK335R','NRCA4_FULL_MASKSWB','NRCA5_FULL_MASKLWB']
grism = ['NRCA5_GRISM_F322W2']
apertures = apertures + coron + grism

#FOR TESTING
apertures = ['NRCA1_FULL','NRCA2_FULL','NRCA3_FULL','NRCA4_FULL','NRCA5_FULL',
             'NRCB1_FULL','NRCB2_FULL','NRCB3_FULL','NRCB4_FULL','NRCB5_FULL']

        
for entry in apertures:
    det = entry[0:5]
    aper_minus_det = entry[6:]

    exp_type = 'NRC_IMAGE'
    pupil = 'ANY'
    if 'MASK' in entry:
        exp_type = 'NRC_CORON'
        m = entry.find('MASK')
        mtype = entry[m+4:]
        if 'R' in mtype:
            pupil = 'MASKRND'
        elif 'B' in mtype:
            pupil = 'MASKBAR'
        else:
            print('Cannot determine type of mask for {}'.format(entry))
            sys.exit()
    elif 'GRISM' in entry:
        exp_type = 'NRC_GRISM'
        pupil = 'GRISMR'

    outname = det + '_' + aper_minus_det + '_distortion.asdf'
    print('running {}'.format(entry))
    ref.create_nircam_distortion(nrc_siaf,det,aper_minus_det,exp_type,pupil,outname)
    

#for entry in siaf:
#    apername = entry['col4']
#    opgs = str(entry['col18'])
#    if apername != 'na' and opgs != 'na' and 'GRISMR' not in apername and 'GRISMC' not in apername:
#        aper_minus_det = apername[6:]
#        det = apername[0:5]
#        try:
#            detnum = np.int(det[-1])
#            outname = 'reffiles_27Oct2016/' + det + '_' + aper_minus_det + '_distortion.asdf'
#            print('running {},{}'.format(apername,opgs))
#            ref.create_nircam_distortion('NIRCam_SIAF_2016-09-29.csv',det,aper_minus_det,opgs,outname)
#        except:
#            #print('skipping {}.'.format(entry))
#            pass
        
