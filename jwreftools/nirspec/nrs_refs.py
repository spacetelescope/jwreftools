import os.path
import argparse
#from asdf import AsdfFile
#from astropy.modeling.models import Mapping, Identity
#from astropy.modeling import models
#from .utils import (linear_from_pcf_det2sky, common_reference_file_keywords,
#                    coeffs_from_pcf)

reftypes = ['camera', 'collimator', 'fpa', 'disperser', 'fore', 'ifuslicer', 'ifupost', 'msa', 'ote', 'wavelenghtrange']


def create_all(model_path):
    pass
    
reftype_mapping = {'camera': camera2asdf, 
                   'collimator': collimator2asdf, 
                   'fpa': fpa2asdf,
                   'disperser': dispoerser2asdf,
                   'fore': fore2asdf, 
                   'ifuslicer': ifuslicer2asdf,
                   'ifupost': ifupost2asdf,
                   'msa': msa2asdf,
                   'ote': ote2asdf,
                   'wavelenghtrange': wave_range2asdf}

  
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Creates NIRSpec reference files in ASDF format.")
    parser.add_argument("mode", type=str, help="'all' or one of reftypes")
    parser.add_argument("model_path", type=str, help="Directory path to NIRSpec Model.")
    parser.add_argument("reffile_path", type=str, help="File path to NIRSpec Reference file.")
    parser.add_argument("output_name", type=str, help="Output file name.")
    res = parser.parse_args()
    if res.mode == "all":
        create_all(res.model_path)
    else:
        if res.reffile_path is None:
            raise ValueError("expect a path to a reference file.")
        reftype_mapping[res.mode](res.reffile_path, res.output_name)

        
   
