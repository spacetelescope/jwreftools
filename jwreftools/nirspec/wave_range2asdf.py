from asdf import AsdfFile
from .utils import common_reference_file_keywords


def wavelength_range(spectral_conf, outname, ref_kw):
    """
    Parameters
    ----------
    spectral_conf : str
        reference file: spectralconfigurations.txt
    outname : str
        output file name
    """
    with open(spectral_conf) as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines][13 :]
    lines = [l.split() for l in lines]
    tree = ref_kw.copy()
    filter_grating = {}
    for l in lines:
        f_g = l[0] + '_' + l[1]
        filter_grating[f_g] = {'order': int(l[2]), 'range': [float(l[3]), float(l[4])]}
    tree['filter_grating'] = filter_grating
    # values in lamp_grating come from private communication with the INS team
    #lamp_grating = {}
    filter_grating['FLAT1_G140M'] = {'order': -1, 'range': [1e-6, 1.8e-6]}
    filter_grating['LINE1_G140M'] = {'order': -1, 'range': [1e-6, 1.8e-6]}
    filter_grating['FLAT1_G140H'] = {'order': -1, 'range': [1e-6, 1.8e-6]}
    filter_grating['LINE1_G140H'] = {'order': -1, 'range': [1e-6, 1.8e-6]}
    filter_grating['FLAT2_G235M'] = {'order': -1, 'range': [1.7e-6, 3.1e-6]}
    filter_grating['LINE2_G235M'] = {'order': -1, 'range': [1.7e-6, 3.1e-6]}
    filter_grating['FLAT2_G235H'] = {'order': -1, 'range': [1.7e-6, 3.1e-6]}
    filter_grating['LINE2_G235H'] = {'order': -1, 'range': [1.7e-6, 3.1e-6]}
    filter_grating['FLAT3_G395M'] = {'order': -1, 'range': [2.9e-6, 5.3e-6]}
    filter_grating['LINE3_G395M'] = {'order': -1, 'range': [2.9e-6, 5.3e-6]}
    filter_grating['FLAT3_G395H'] = {'order': -1, 'range': [2.9e-6, 5.3e-6]}
    filter_grating['LINE3_G395H'] = {'order': -1, 'range': [2.9e-6, 5.3e-6]}
    filter_grating['REF_G140M'] = {'order': -1, 'range': [1.3e-6, 1.7e-6]}
    filter_grating['REF_G140H'] = {'order': -1, 'range': [1.3e-6, 1.7e-6]}
    filter_grating['TEST_MIRROR'] = {'order': -1, 'range': [0.6e-6, 5.3e-6]}
    #for grating in ["G140H", "G140M", "G235H", "G235M", "G395H", "G395M", "MIRROR"]:
        #lamp_grating['FLAT4_{0}'.format(grating)] = {'order': -1, 'range': [0.7e-6, 1.2e-6]}
    #for grating in ["G140H", "G140M", "G235H", "G235M", "G395H", "G395M", "MIRROR"]:
        #lamp_grating['LINE4_{0}'.format(grating)] = {'order': -1, 'range': [0.6e-6, 5.3e-6]}
    #tree['lamp_grating'] = lamp_grating
    fasdf = AsdfFile()

    fasdf.tree = tree
    fasdf.add_history_entry("Build 6")
    fasdf.write_to(outname)
    return fasdf

if __name__ == '__main__':
    import argpars
    parser = argpars.ArgumentParser(description="Creates NIRSpec 'wavelengthrange' reference file in ASDF format.")
    parser.add_argument("wave_range_file", type=str, help="Spectral configurations file.")
    parser.add_argument("output_name", type=str, help="Output file name")
    res = parser.parse_args()
    if res.output_name is None:
        output_name = "nirspec_wavelength.asdf"
    else:
        output_name = res.output_name

    ref_kw = common_reference_file_keywords("WAVELENGTHRANGE", "NIRSPEC Spectral Configurations - CDP4")
    wavelength_range(wave_range_file, output_name, ref_kw)
