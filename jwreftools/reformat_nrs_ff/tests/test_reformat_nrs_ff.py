"""
Test for reformat_nrs_ff
"""
import os
import pytest
from tempfile import TemporaryDirectory

import numpy as np

from astropy.io import fits

from .. import reformat_nrs_ff


@pytest.fixture
def data_file():
    # Create a dummy SFLAT file for NRS_FIXEDSLIT data.
    ifd = fits.HDUList(fits.PrimaryHDU())
    ifd[0].header['telescop'] = 'JWST'
    ifd[0].header['detector'] = 'NRS2'
    ifd[0].header['insrtume'] = 'NIRSPEC'   # keyword intentionally misspelled
    ifd[0].header['filter'] = 'OPAQUE'
    ifd[0].header['detector'] = 'NRS2'
    ifd[0].header['grating'] = 'G235H'
    ifd[0].header['useafter'] = '2016-06-28'
    ifd[0].header['reftype'] = 'CAAFLAT2'
    ifd[0].header['descrip'] = 'Flat-field reference file for Fixed Slits'
    ifd[0].header['pedigree'] = 'GROUND'

    sci_data = np.arange(5 * 5, dtype=np.float32).reshape((5, 5)) + 2
    err_data = np.arange(5 * 5, dtype=np.float32).reshape((5, 5)) + 1
    dq_data = np.arange(5 * 5, dtype=np.float32).reshape((5, 5))
    ifd.append(fits.ImageHDU(data=sci_data, name='SCI'))
    ifd.append(fits.ImageHDU(data=err_data, name='ERR'))
    ifd.append(fits.ImageHDU(data=dq_data, name='DQ'))
    del sci_data, err_data, dq_data

    # DQ_DEF
    col = []
    col.append(fits.Column(name='BIT', format='1J', array=[0, 1, 2]))
    col.append(fits.Column(name='VAL', format='1J', array=[1, 2, 4]))
    col.append(fits.Column(name='NAME', format='20A',
                           array=['DO_NOT_USE',
                                  'UNRELIABLE_FLAT',
                                  'NO_FLAT_FIELD']))
    col.append(fits.Column(name='DESCRIPTION', format='40A',
                           array=['Bad pixel. Do not use.',
                                  'Flat variance large',
                                  'Flat field cannot be measured']))
    dq_def = fits.BinTableHDU.from_columns(fits.ColDefs(col))
    dq_def.header['extname'] = 'DQ_DEF'
    ifd.append(dq_def)

    # Fast-variation tables.
    # This wavelength column will be used for all five tables.
    wl_col = fits.Column(name='WAVELENGTH', unit='MICROMETER',
                         format='1E',
                         array=[1.66, 2.16332, 2.66665, 3.16997])
    col = []
    col.append(wl_col)
    col.append(fits.Column(name='DATA', unit='UNITLESS',
                           format='1E',
                           array=[0.26492, 0.26492, 1.63155, 1.63155]))
    slit = fits.BinTableHDU.from_columns(fits.ColDefs(col))
    slit.header['extname'] = 'SLIT_A_200_1'
    ifd.append(slit)
    col = []
    col.append(wl_col)
    col.append(fits.Column(name='DATA', unit='UNITLESS',
                           format='1E',
                           array=[0.27124, 0.27124, 1.64399, 1.64399]))
    slit = fits.BinTableHDU.from_columns(fits.ColDefs(col))
    slit.header['extname'] = 'SLIT_A_200_2'
    ifd.append(slit)
    col = []
    col.append(wl_col)
    col.append(fits.Column(name='DATA', unit='UNITLESS',
                           format='1E',
                           array=[0.52911, 0.52911, 3.54949, 3.54949]))
    slit = fits.BinTableHDU.from_columns(fits.ColDefs(col))
    slit.header['extname'] = 'SLIT_A_400'
    ifd.append(slit)
    col = []
    col.append(wl_col)
    col.append(fits.Column(name='DATA', unit='UNITLESS',
                           format='1E',
                           array=[0.42394, 0.42394, 1.64300, 1.64300]))
    slit = fits.BinTableHDU.from_columns(fits.ColDefs(col))
    slit.header['extname'] = 'SLIT_B_200'
    ifd.append(slit)
    col = []
    col.append(wl_col)
    col.append(fits.Column(name='DATA', unit='UNITLESS',
                           format='1E',
                           array=[2.29804, 2.29804, 15.0267, 15.0267]))
    slit = fits.BinTableHDU.from_columns(fits.ColDefs(col))
    slit.header['extname'] = 'SLIT_A_1600'
    ifd.append(slit)

    """The exposure type (e.g. NRS_FIXEDSLIT) and optical path part
       (FFLAT, SFLAT, or DFLAT) are inferred from the file name.
       For SFLAT (only), if keyword FILTER is "OPAQUE", the relevant
       filter name for science data will be determined from the
       "FLATi" (i = 1..5) string in the file name.
       nirspec_FS_sflat_G235H_OPAQUE_FLAT2_nrs2_f_01.01.fits
    """

    with TemporaryDirectory() as dirname:
        file_path = os.path.join(
                dirname,
                "nirspec_FS_sflat_G235H_OPAQUE_FLAT2_nrs2_f_01.01.fits")
        ifd.writeto(file_path)
        ifd.close()
        yield file_path


def test_reformat_nrs_ff(data_file):
    # This function reformats the IDT version to be JWST pipeline compatible.

    words = os.path.split(data_file)
    prefix = os.path.join(words[0], "xyz_")
    output_file = prefix + words[1]

    reformat_nrs_ff.process_files(data_file, prefix=prefix, verbose=False)

    ofd = fits.open(output_file)
    phdr = ofd[0].header
    assert phdr['telescop'] == 'JWST'
    assert phdr['detector'] == 'NRS2'
    assert phdr['instrume'] == 'NIRSPEC'
    assert phdr['useafter'] == '2016-06-28T00:00:00'
    assert phdr['exp_type'] == 'NRS_FIXEDSLIT'
    assert phdr['reftype'] == 'SFLAT'
    assert phdr['filter'] == 'F170LP'

    sci_ref = np.array([[26., 21., 16., 11., 6.],
                        [25., 20., 15., 10., 5.],
                        [24., 19., 14.,  9., 4.],
                        [23., 18., 13.,  8., 3.],
                        [22., 17., 12.,  7., 2.]], dtype=np.float32)
    assert np.allclose(ofd[("sci", 1)].data, sci_ref, rtol=0., atol=1.e-6)
    del sci_ref

    dq_ref = np.array([[24, 19, 14, 9, 4],
                       [23, 18, 13, 8, 3],
                       [22, 17, 12, 7, 2],
                       [21, 16, 11, 6, 1],
                       [20, 15, 10, 5, 0]], dtype=np.int32)
    assert np.all(np.equal(ofd[("dq", 1)].data, dq_ref))
    del dq_ref

    err_ref = np.array([[25., 20., 15., 10., 5.],
                        [24., 19., 14., 9., 4.],
                        [23., 18., 13., 8., 3.],
                        [22., 17., 12., 7., 2.],
                        [21., 16., 11., 6., 1.]], dtype=np.float32)
    assert np.allclose(ofd[("err", 1)].data, err_ref, rtol=0., atol=1.e-6)
    del err_ref

    slit_name_ref = np.array(['S200A1', 'S200A2', 'S400A1', 'S200B1',
                              'S1600A1'])
    nelem_ref = np.array([4, 4, 4, 4, 4], dtype=np.int32)
    wavelength_ref = np.array([[1.66, 2.16332, 2.66665, 3.16997],
                               [1.66, 2.16332, 2.66665, 3.16997],
                               [1.66, 2.16332, 2.66665, 3.16997],
                               [1.66, 2.16332, 2.66665, 3.16997],
                               [1.66, 2.16332, 2.66665, 3.16997]],
                              dtype=np.float32)
    data_ref = np.array([[ 0.26492, 0.26492, 1.63155, 1.63155],
                         [ 0.27124, 0.27124, 1.64399, 1.64399],
                         [ 0.52911, 0.52911, 3.54949, 3.54949],
                         [ 0.42394, 0.42394, 1.643,   1.643  ],
                         [ 2.29804, 2.29804, 15.0267, 15.0267 ]],
                        dtype=np.float32)

    fast_var = ofd[("fast_variation", 1)].data
    slit_name = fast_var.field("slit_name")
    nelem = fast_var.field("nelem")
    wavelength = fast_var.field("wavelength")
    data = fast_var.field("data")

    for i in range(len(slit_name_ref)):
        assert slit_name[i] == slit_name_ref[i]

    assert np.all(np.equal(nelem, nelem_ref))

    assert np.allclose(wavelength, wavelength_ref, rtol=0., atol=1.e-6)

    assert np.allclose(data, data_ref, rtol=0., atol=1.e-6)

    ofd.close()
