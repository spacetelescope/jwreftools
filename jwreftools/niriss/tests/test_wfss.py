import glob
import pytest
import numpy as np
from numpy.testing import assert_allclose
from astropy.modeling.models import *
from astropy.modeling.models import math as astmath
from gwcs.wcstools import grid_from_bounding_box

from jwst import datamodels
from jwst.transforms import models

import grismconf
from .. import niriss_grism_reffiles as ngr

from . import data


def create_specwcs(conf, outname='test.asdf'):
    filter_name, pupil = _get_configuration(conf)
    ngr.create_grism_config(conf,
               filter_name,
               pupil,
               354.222229,
               author="STScI",
               history="NIRISS Grism Parameters",
               outname=outname)
    print("Created {outname} from ", conf)


def get_conf_files():
    data_path = data.__path__[0]
    return glob.glob(data_path + '/*.conf')


def _get_configuration(conf):
    name = conf.split("/")[-1]
    _, pupil, filter_name = name.split('_')
    pupil = pupil.split('.')[0].upper()
    return filter_name, pupil

#pytest.mark.parametrize('order', ['-1', '0', '+1', '2', '3'])
def test_spec_vs_grismconf():
    # "/Users/dencheva/dev/GRISM_NIRISS/V3/NIRISS_F200W_GR150r.conf"
    conf_files = get_conf_files()
    for i, conf in enumerate(conf_files):
        print("working on file ", conf)
        #filter_name, pupil = _get_configuration(conf)
        name = f"jwst_niriss_specwcs_{i}.asdf"
        create_specwcs(conf, name)
        C = grismconf.Config(conf)
        x0 = 500.5
        y0 = 600.1
        order='+1'

        t = np.linspace(0, 1, 10)
        wavelength = C.DISPL(order,x0,y0,t)
        grismconf_dx = C.DISPX(order,x0,y0,t)
        grismconf_dy = C.DISPY(order,x0,y0,t)

        tspec = datamodels.open(name)
        #tspec.fwcpos_ref = 354.222229
        orders = tspec.orders
        ind = orders.index(int(order))
        xmodels=tspec.dispx[ind]
        ymodels=tspec.dispy[ind]
        pipe_dx = xmodels[0](x0, y0) + t * xmodels[1](x0, y0) + t**2 * xmodels[2](x0, y0)
        pipe_dy = ymodels[0](x0, y0) + t * ymodels[1](x0, y0) + t**2 * ymodels[2](x0, y0)

        assert_allclose(pipe_dx, grismconf_dx)
        assert_allclose(pipe_dy, grismconf_dy)

        # Load the Grism Configuration file
        # C = grismconf.Config("NIRISS_F090W_GR150C.conf")
        # Compute the t values corresponding to the exact offsets
        ts = C.INVDISPX(order, x0, y0, grismconf_dx, t)
        assert_allclose(t, ts)

        # Compute the dys values for the same pixels
        dys = C.DISPY("+1", x0, y0, ts)
        so = np.argsort(pipe_dx)
        tab = Tabular1D(pipe_dx[so], t[so], bounds_error=False, fill_value=None)

        pipe_ts = tab(pipe_dx)
        assert_allclose(pipe_ts, ts)
        #ts - pipe_ts

        dyn = ymodels[0](x0, y0) + pipe_ts * ymodels[1](x0, y0) + pipe_ts**2 * ymodels[2](x0, y0)
        assert_allclose(dys, dyn)

        lmodels = tspec.displ[ind]
        lam = lmodels(pipe_ts)

        # Compute wavelength of each of the pixels
        wavs = C.DISPL("+1", x0, y0, ts)
        assert_allclose(lam, wavs)
        # print("\n dys", dys)
        # print("\n wavs", wavs)
