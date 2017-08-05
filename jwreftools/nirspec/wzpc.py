from astropy.io import fits
from astropy.modeling import models
from astropy import wcs
import numpy as np
from jwst.datamodels import WZPCModel


def _wzpc2asdf(wzpcfile, author, description, useafter):

    f = fits.open(wzpcfile)
    width = f[1].header['width']
    name = f[0].header['COMPNAME']
    data = f[1].data
    var = f[2].data
    w = wcs.WCS(f[1].header)
    f.close()
    
    y, x = np.mgrid[: data.shape[0], : data.shape[1]]
    X, Y = w.all_pix2world(x, y, 1)
    tab = models.Tabular2D(points=(X[0], Y[:,0]), lookup_table=data)
    aps=[{'name': name, 'variance': var, 'zero_point_offset': tab, 'width': width}]
