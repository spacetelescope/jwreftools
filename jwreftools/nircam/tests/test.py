import os.path
import numpy as np
from numpy.testing import assert_allclose

import pytest
import grismconf

from astropy.modeling.models import Tabular1D
from jwreftools.nircam import nircam_grism_reffiles as nrc
from jwst import datamodels
from . import data

data_path = data.__path__[0]

pytest.mark.parametrize('order', ['+1', '+2'])
def test_vs_grismconf(order):
    x0 = 500.5
    y0 = 600.1
    order = order

    tspec = datamodels.open(os.path.join(data_path, 'test_nircam_specwcs.asdf'))
    conf = os.path.join(data_path, 'NIRCAM_F360M_modA_C.conf')
    C = grismconf.Config(conf)

    t = np.linspace(0, 1, 10)
    wavelength = C.DISPL(order,x0,y0,t)
    grismconf_dx = C.DISPX(order,x0,y0,t)
    grismconf_dy = C.DISPY(order,x0,y0,t)
    orders = tspec.orders
    ind = orders.index(int(order))
    xmodels=tspec.dispx[ind]
    ymodels=tspec.dispy[ind]
    if xmodels[0].n_inputs == 2:
        pipe_dx = xmodels[0](x0, y0) + t * xmodels[1](x0, y0) + t**2 * xmodels[2](x0, y0)
        pipe_dy = ymodels[0](x0, y0) + t * ymodels[1](x0, y0) + t**2 * ymodels[2](x0, y0)
    elif xmodels[0].n_inputs == 1:
        pipe_dx = xmodels[0](t)
        pipe_dy = ymodels[0](t)

    assert_allclose(pipe_dx, grismconf_dx)
    assert_allclose(pipe_dy, grismconf_dy)

    # Load the Grism Configuration file
    # C = grismconf.Config("NIRISS_F090W_GR150C.conf")
    # Compute the t values corresponding to the exact offsets
    ts = C.INVDISPX(order, x0, y0, grismconf_dx, t)
    assert_allclose(t, ts)

    dys = C.DISPY(order, x0, y0, ts)
    so = np.argsort(pipe_dx)
    tab = Tabular1D(pipe_dx[so], t[so], bounds_error=False, fill_value=None)

    pipe_ts = tab(pipe_dx)
    assert_allclose(pipe_ts, ts)

    if ymodels[0].n_inputs == 2:
        dyn = ymodels[0](x0, y0) + pipe_ts * ymodels[1](x0, y0) + pipe_ts**2 * ymodels[2](x0, y0)

    elif ymodels[0].n_inputs == 1:
        dyn = ymodels[0](pipe_ts)
    assert_allclose(dys, dyn)

    lmodels = tspec.displ[ind]
    if lmodels[0].n_inputs == 2:
        lam = lmodels[0](x0, y0) + pipe_ts * lmodels[1](x0, y0) + pipe_ts**2 * lmodels[2](x0, y0)
    elif lmodels[0].n_inputs == 1:
        lam = lmodels[0](pipe_ts)

    # Compute wavelength of each of the pixels
    wavs = C.DISPL(order, x0, y0, ts)
    assert_allclose(lam, wavs)
