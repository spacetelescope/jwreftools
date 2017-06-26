import nircam_reftools as ref
from astropy.io import ascii
import numpy as np

siaf_file = 'NIRCam_SIAF-OSS_Summary_20170329_MMA.csv'
nrc_siaf = 'NIRCam_SIAF_2017-03-28.csv'

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

#FOR TESTING
grism_apertures = ['NRCA5_GRISM_F322W2']


#necessary IMAGING metadata-----------------------------------------------
imaging_apertures = ['NRCA1_FULL','NRCA2_FULL','NRCA3_FULL','NRCA4_FULL'
                     ,'NRCA5_FULL',
                     'NRCB1_FULL','NRCB2_FULL','NRCB3_FULL','NRCB4_FULL'
                     ,'NRCB5_FULL']
sw_imaging_pupil = ['CLEAR','F162M','F164N','GDHS0','GDHS60','WLM8'
                    ,'WLP8','PINHOLES','MASKIPR']
lw_imaging_pupil = ['CLEAR','F323N','F405N','F466N','F470N','PINHOLES'
                    ,'MASKIPR']

sw_imaging_subarr = ['FULL','SUB64P','SUB160','SUB160P','SUB320','SUB400P'
                     ,'SUB640','SUB8FP1A','SUB8FP1B','SUB64FP1A','SUB64FP1B'
                     ,'SUB96DHSPILA','SUB96DHSPILB']
lw_imaging_subarr = ['FULL','SUB64P','SUB160','SUB160P','SUB320','SUB400P'
                     ,'SUB640','SUB32TATS','SUB32TATSGRISM']

imaging_exptype = ['NRC_IMAGE','NRC_TSIMAGE','NRC_FLAT','NRC_LED'
                   ,'NRC_WFSC']
#-------------------------------------------------------------------------

#necessary CORONAGRAPHY metadata------------------------------------------

#SW ROUND, A2
sw_around_apertures = 'NRCA2_FULL_MASK210R'
sw_around_pupil = ['MASKRND']
sw_around_subarr = ['FULL','SUB640A210R','SUBNDA210R','SUBFSA210R']

#SW BAR, A4
sw_abar_apertures = 'NRCA4_FULL_MASKSWB'
sw_abar_pupil = ['MASKBAR']
sw_abar_subarr = ['FULL','SUB640ASWB','SUBFSASWB','SUBNDASWBL','SUBNDASWBS']

#LW BAR, A5
lw_abar_apertures = 'NRCA5_FULL_MASKLWB'
lw_abar_pupil = ['MASKBAR']
lw_abar_subarr = ['FULL','SUB320ALWB','SUBFSALWB','SUBNDALWBL','SUBNDALWBS']

#LW ROUND, A5
lw_around_apertures = 'NRCA5_FULL_MASK335R'
lw_around_pupil = ['MASKRND']
lw_around_subarr = ['FULL','SUB320A335R','SUB320A430R','SUBFSA335R'
                    ,'SUBFSA430R','SUBNDA335R','SUBNDA430R']

#SW BAR, B3
#sw_bbar_apertures = ['NRCB3_MASKSWB']
#sw_bbar_pupil = ['MASKSWBAR']
#sw_bbar_subarr = ['FULL','SUB640BSWB','SUBNDBSWBL','SUBNDBSWBS']

#SW ROUND, B1
#sw_bround_apertures = ['NRCB1_MASK210R']
#sw_bround_pupil = ['MASKSWRND']
#sw_bround_subarr = ['FULL','SUB640B210R','SUBNDB210R']

#LW BAR, B5
#lw_bbar_apertures = ['NRCB5_MASKLWB']
#lw_bbar_pupil = ['MASKLWB']
#lw_bbar_subarr = ['FULL','SUBNDBLWBL','SUBNDBLWBS']

#LW ROUND, B5
#lw_bround_apertures = ['NRCB5_MASK335R']
#lw_bround_pupil = ['MASKLWRND']
#lw_bround_subarr = ['FULL','SUB320B335R','SUB320B430R','SUBNDB335R','SUBNDB430R']

coron_exptypes = ['NRC_CORON','NRC_TACONFIRM','NRC_TACQ']


#coron_apertures = ['NRCA2_MASK210R','NRCA5_MASK335R',
#                   'NRCA4_MASKSWB','NRCA5_MASKLWB',
#                   'NRCB1_MASK210R','NRCB5_MASK335R',
#                   'NRCB3_MASKSWB','NRCB5_MASKLWB']
#swbar_pupil = ['MASKSWBAR']
#swbar_subarr = ['FULL','SUB640ASWB','SUB640BSWB','SUBFSASWB'
#                ,'SUBNDASWBL','SUBNDASWBS','SUBNDBSWBL','SUBNDBSWBS']
#swbar_exptype = ['NRC_CORON','NRC_TACONFIRM','NRC_TACQ']

#swrnd_pupil = ['MASKSWRND']
#swrnd_subarr = ['FULL','SUB640A210R','SUB640B210R','SUBFSA210R'
#                ,'SUBNDA210R','SUBNDB210R']
#swrnd_exptype = ['NRC_CORON','NRC_TACONFIRM','NRC_TACQ']

#lwbar_pupil = ['MASKLWBAR']
#lwbar_subarr = ['FULL','SUB320ALWB','SUB320BLWB','SUBFSALWB','SUBNDALWBL'
#                ,'SUBNDALWBS','SUBNDBLWBL','SUBNDBLWBS']
#lwbar_exptype = ['NRC_CORON','NRC_TACONFIRM','NRC_TACQ']

#lwrnd_pupil = ['MASKLWRND']
#lwrnd_subarr = ['FULL','SUB320A335R','SUB320A430R','SUB320B335R'
#                ,'SUB320B430R','SUBFSA335R','SUBFSA430R','SUBNDA335R'
#                ,'SUBNDA430R','SUBNDB335R','SUBNDB430R']
#lwrnd_exptype = ['NRC_CORON','NRC_TACONFIRM','NRC_TACQ']
#-------------------------------------------------------------------------


#RUN IMAGING
for entry in imaging_apertures:
    det = entry[0:5]
    aper_minus_det = entry[6:]

    if np.int(det[-1]) < 5:
        pupil = sw_imaging_pupil
        subarr = sw_imaging_subarr
    else:
        pupil = lw_imaging_pupil
        subarr = lw_imaging_subarr
        
    exp_type = imaging_exptype

    outname = det + '_' + aper_minus_det + '_distortion.asdf'
    #print('running {}'.format(entry))
    #print('det is {}'.format(np.int(det[-1])))
    #print('pupil is {}'.format(pupil))
    ref.create_nircam_distortion(nrc_siaf,det,aper_minus_det,outname,pupil,subarr,exp_type)
    

#RUN CORONAGRAPHY - NOTE!!! CORONAGRAPHY IS ONLY SUPPORTED IN MODULE A!

c_apertures = [sw_around_apertures,sw_abar_apertures,lw_around_apertures,lw_abar_apertures]
c_pupil = [sw_around_pupil,sw_abar_pupil,lw_around_pupil,lw_abar_pupil]
c_subarr = [sw_around_subarr,sw_abar_subarr,lw_around_subarr,lw_abar_subarr]

for ap,pup,sub in zip(c_apertures,c_pupil,c_subarr):
    det = ap[0:5]
    aper_minus_det = ap[6:]
    outname = ap + '_distortion.asdf'
    ref.create_nircam_distortion(nrc_siaf,det,aper_minus_det,outname,pup,sub,coron_exptypes)
        

    
