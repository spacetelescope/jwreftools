from astropy.io import fits
from astropy.modeling import models
from astropy import wcs
import numpy as np
from jwst.datamodels import WZPCModel


ap_names_map = {'A200_1': 'S200A1',
                'A200_2': 'S200A2',
                'A400': 'S400A1',
                'A1600': 'S1600A1',
                'B200': 'S200B1'
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
    X, Y = w.all_pix2world(x, y, 1)
    tab = models.Tabular2D(points=(X[0], Y[:,0]), lookup_table=data)
    aperture = {'aperture_name': name, 'variance': var, 'zero_point_offset': tab, 'width': width}
    return aperture


def create_wzpc_refs(wzpc_files, outname=None, author=None, description=None, useafter="2015-11-01"):
    """
    Create WZPC reference files.

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
    model = WZPCModel()
    aps = []
    if isinstance(wzpc_files, list):
        # Create a reference file for the Fixed Slits mode.
        f0 = wzpc_files[0]
        model.meta.exposure.type = "NRS_FIXEDSLIT"
        for f in wzpc_files:
            aps.append(_wzpc2asdf(f, author=author, description=description, useafter=useafter))
    elif isinstance(wzpc_files, str):
        model.meta.exposure.type = "NRS_MSASPEC"
        f0 = wzpc_files
        aps.append(_wzpc2asdf(wzpc_files, author=author, description=description, useafter=useafter))
    else:
        raise ValueError("Invalid input - expected a string or a list of strings.")

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

    if outname is None:
        outname = "nirspec_wzpc.asdf"
    model.to_asdf(outname)
