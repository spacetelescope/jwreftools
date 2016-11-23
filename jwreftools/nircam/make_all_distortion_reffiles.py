import nircam_reftools as ref
from astropy.io import ascii
import numpy as np

siaf_file = 'NIRCam_SIAF-OSS_Summary_20161117_MMA.csv'

siaf = ascii.read(siaf_file,data_start=18)
siaf['col4'].fill_value = 'na'
siaf['col18'].fill_value = 'na'
siaf = siaf.filled()

for entry in siaf:
    apername = entry['col4']
    opgs = str(entry['col18'])
    if apername != 'na' and opgs != 'na' and 'GRISMR' not in apername and 'GRISMC' not in apername:
        aper_minus_det = apername[6:]
        det = apername[0:5]
        try:
            detnum = np.int(det[-1])
            outname = 'reffiles_27Oct2016/' + det + '_' + aper_minus_det + '_distortion.asdf'
            print('running {},{}'.format(apername,opgs))
            ref.create_nircam_distortion('NIRCam_SIAF_2016-09-29.csv',det,aper_minus_det,opgs,outname)
        except:
            #print('skipping {}.'.format(entry))
            pass
        
