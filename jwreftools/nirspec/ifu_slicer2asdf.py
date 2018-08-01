import datetime
import os.path
import numpy as np
from astropy.modeling import models
from astropy.io import fits
from jwst.datamodels import IFUSlicerModel
from asdf.tags.core import Software, HistoryEntry


__all__ = ["create_ifuslicer_reference", "ifu_slicer2asdf"]


def ifu_slicer2asdf(ifuslicer, author, description, useafter):
    """
    Create an asdf reference file with the MSA description.

    ifu_slicer2asdf("IFU_slicer.sgd", "ifu_slicer.asdf")

    Parameters
    ----------
    ifuslicer : str
        A fits file with the IFU slicer description
    outname : str
        Name of output ASDF file.
    """
    f = fits.open(ifuslicer)
    data = f[1].data
    header = f[1].header
    shiftx = models.Shift(header['XREF'], name='ifuslicer_x')
    shifty = models.Shift(header['YREF'], name='ifuslicer_y')
    rot = models.Rotation2D(np.rad2deg(header['ROT']), name='ifuslicer_rotate')
    model = rot | shiftx & shifty
    f.close()

    slicer_model = IFUSlicerModel()
    slicer_model.model = model
    slicer_model.data = data
    slicer_model.meta.author = author
    slicer_model.meta.description = description
    slicer_model.meta.useafter = useafter
    slicer_model.meta.pedigree = "GROUND"
    return slicer_model


def create_ifuslicer_reference(ifuslicer_refname, output_name, author=None, description=None, useafter=None):
    with fits.open(ifuslicer_refname) as f:
        auth = f[0].header['AUTHOR']
        descrip = f[0].header['DESCR']
        date = f[0].header['DATE']
    if author is None:
        author = auth
    if description is None:
        description = descrip
    if useafter is None:
        useafter = date

    try:
        model = ifu_slicer2asdf(ifuslicer_refname, author, description, useafter)
    except:
        raise Exception("IFUSLICER file was not created.")
    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    model.history.append(entry)
    model.to_asdf(output_name)
    model.validate()
