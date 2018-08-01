import datetime
import numpy as np
from jwst.datamodels import WavelengthrangeModel
from astropy import units as u
from asdf.tags.core import Software, HistoryEntry


__all__ = ["create_wavelengthrange_reference", "wavelength_range"]


def wavelength_range(spectral_conf, author, description, useafter):
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
    lines = [l.strip() for l in lines]
    for index, line in enumerate(lines):
        if 'CONFIGURATIONS' in line:
            break
    lines = np.array([l.split() for l in lines[index + 2 :]]).T
    filt, grat, order, wmin, wmax, samping = lines
    filter_grating = [f+'_'+g for f, g in zip(filt, grat)]
    order = [int(i) for i in order]
    wave_range = [[mini, maxi] for mini, maxi in zip(
        wmin.astype(np.float), wmax.astype(np.float))]

    # in addition
    addon_list = ['FLAT1_G140H', 'FLAT2_G235H', 'LINE3_G395H', 'REF_G140H',
                  'REF_G140M', 'FLAT2_G235M', 'LINE1_G140H', 'FLAT3_G395M',
                  'LINE1_G140M', 'FLAT3_G395H', 'FLAT1_G140M', 'LINE2_G235M',
                  'LINE3_G395M', 'LINE2_G235H', 'TEST_MIRROR']

    filter_grating.extend(addon_list)
    order.extend([-1] * len(addon_list))
    wave_range.extend([[1e-06, 1.8e-06], [1.7e-06, 3.1e-06], [2.9e-06, 5.3e-06], [1.3e-06, 1.7e-06],
                       [1.3e-06, 1.7e-06], [1.7e-06, 3.1e-06], [1e-06, 1.8e-06], [2.9e-06, 5.3e-06],
                       [1e-06, 1.8e-06], [2.9e-06, 5.3e-06], [1e-06, 1.8e-06], [1.7e-06, 3.1e-06],
                       [2.9e-06, 5.3e-06], [1.7e-06, 3.1e-06], [6e-07, 5.3e-06]])

    wr_model = WavelengthrangeModel()
    wr_model.waverange_selector = filter_grating
    wr_model.wavelengthrange = wave_range
    wr_model.order = order
    wr_model.meta.wavelength_units = u.m
    wr_model.meta.author = author
    wr_model.meta.description = description
    wr_model.meta.useafter = useafter
    wr_model.meta.pedigree = "GROUND"
    wr_model.meta.instrument.name = "NIRSPEC"
    wr_model.meta.instrument.p_detector = "NRS1|NRS2|"
    wr_model.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
    NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
    NRS_MIMF|NRS_BOTA|NRS_LAMP|NRS_BRIGHTOBJ|"
    wr_model.meta.exposure.type = "N/A"
    return wr_model

def create_wavelengthrange_reference(wave_range_file, output_name, author=None, description=None, useafter=None):
    f = open(wave_range_file)
    lines = [l.strip() for l in f.readlines()]
    f.close()
    for i, line in enumerate(lines):
        if 'AUTHOR' in line:
            auth = lines[i + 1]
            continue
        elif 'DESCRIPTION' in line:
            descrip = lines[i + 1]
            continue
        elif 'DATE' in line:
            date = lines[i + 1]
            continue

    if author is None:
        author = auth
    if description is None:
        description = descrip
    if useafter is None:
        useafter = date

    wr_model = wavelength_range(wave_range_file, author, description, useafter)
    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    wr_model.history.append(entry)
    wr_model.to_asdf(output_name)
    wr_model.validate()
