"""
Create a Nirspec WAVECORR reference file.

Two reference files should be delivered - one for the MSA
and one for the fixed slits. The IDT team delivers a separate file for each slit.
One CRDS ref file is made using the five IDT files.

Examples
--------

Create the fixed slit ref file.

>>> from jwreftools.nirspec import wavecorr2asdf
>>> l = ['jwst-nirspec-slit-a1600.wzprf.fits', 'jwst-nirspec-slit-a200-1.wzprf.fits', 
         'jwst-nirspec-slit-a200-2.wzprf.fits', 'jwst-nirspec-slit-a400.wzprf.fits', 
         'jwst-nirspec-slit-b200.wzprf.fits']
>>> wavecorr2asdf.create_wavecorr_refs(l, author='ESA', outname='nirspec_fs_wavecorr.asdf')

Create the MSA ref file:

>>> wavecorr2asdf.create_wavecorr_refs('jwst-nirspec-mos.wzprf.fits', author='ESA',
        outname='nirspec_mos_wavecorr.asdf')

"""
import datetime
from astropy.io import fits
from astropy.modeling import models
from astropy import wcs
import numpy as np
from asdf.tags.core import HistoryEntry, Software
from jwst.datamodels import WaveCorrModel


ap_names_map = {'A200_1': 'S200A1',
                'A200_2': 'S200A2',
                'A400': 'S400A1',
                'A1600': 'S1600A1',
                'B200': 'S200B1',
                'MOS': 'MOS'
                }


def _wzpc2asdf(wzpcfile, author, description, useafter):

    f = fits.open(wzpcfile)
    width = f[1].header['width']
    name = ap_names_map[f[0].header['COMPNAME']]
    data = f[1].data
    var = f[2].data
    w = wcs.WCS(f[1].header)
    f.close()
    
    y, x = np.mgrid[: data.shape[0], : data.shape[1]]
    # The WCS in the current ref files is 0-based.
    X, Y = w.all_pix2world(x, y, 0)
    tab = models.Tabular2D(points=(X[0], Y[:,0]), lookup_table=data)
    aperture = {'aperture_name': name, 'variance': var, 'zero_point_offset': tab, 'width': width}
    return aperture


def create_wavecorr_refs(wzpc_files, outname=None, author=None, description=None, useafter="2015-11-01"):
    """
    Create WAVECORR reference files (Nirspec wavelength zero-point correction).

    Parameters
    ----------
    wzpc_files : list or str
        A list of *slit-*wzprf.fits file names to create the Fixed Slits reference file.
        A file name to create the MSA reference file.
    author : str
        Author. If None it will be read from the first file.
    description : str
        Description of the file. If None will be read from the header.
    useafter : str
        Useafter date.

    """
    model = WaveCorrModel()
    aps = []
    if isinstance(wzpc_files, list):
        # Create a reference file for the Fixed Slits mode.
        model.meta.exposure.type = "NRS_FIXEDSLIT"
        model.meta.exposure.p_exptype = "NRS_FIXEDSLIT|NRS_BRIGHTOBJ|"
        for f in wzpc_files:
            aps.append(_wzpc2asdf(f, author=author, description=description, useafter=useafter))
        if description is None:
            description = "Wavelength zero-point reference file for Nirspec fixed slits, computed using a simple toy model."
    elif isinstance(wzpc_files, str):
        model.meta.exposure.type = "NRS_MSASPEC"
        wzpc_files = [wzpc_files]
        aps.append(_wzpc2asdf(wzpc_files[0], author=author, description=description, useafter=useafter))
    else:
        raise ValueError("Invalid input - expected a string or a list of strings.")
    f0 = wzpc_files[0]
    model.apertures = aps
    if author is None:
        author = fits.getval(f0, 'AUTHOR')
    if description is None:
        description = fits.getval(f0, 'DESCRIP')
    pedigree = fits.getval(f0, 'PEDIGREE')

    model.meta.author = author
    model.meta.pedigree = pedigree
    model.meta.description = description
    model.meta.useafter = useafter
    model.meta.date = fits.getval(f0, 'date')
    model.meta.origin = fits.getval(f0, 'author')
    model.meta.instrument.p_detector = "NRS1|NRS2|"

    entry = HistoryEntry({'description': "NIRSPEC wavelength zero-point correction.", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    model.history = [entry]
    if outname is None:
        outname = "nirspec_wavecorr.asdf"
    model.to_asdf(outname)
