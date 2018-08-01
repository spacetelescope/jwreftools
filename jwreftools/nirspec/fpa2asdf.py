import datetime
import numpy as np
from jwst.datamodels import FPAModel
from astropy.modeling import models
from asdf.tags.core import HistoryEntry, Software

__all__ = ["create_fpa_reference", "fpa2asdf"]


def fpa2asdf(fpafile, author, description, useafter):
    """
    Create an asdf reference file with the FPA description.

    The CDP2 delivery includes a fits file - "FPA.fpa" which is the
    input to this function. This file is converted to asdf and is a
    reference file of type "FPA".

    fpa2asdf('Ref_Files/CoordTransform/Description/FPA.fpa', 'fpa.asdf')

    Parameters
    ----------
    fpafile : str
        A fits file with FPA description (FPA.fpa)
    outname : str
        Name of output ASDF file.
    """
    with open(fpafile) as f:
        lines = [l.strip() for l in f.readlines()]

    # NRS1
    ind = lines.index("*SCA491_PitchX")
    nrs1_pitchx = float(lines[ind+1])
    ind = lines.index("*SCA491_PitchY")
    nrs1_pitchy = float(lines[ind+1])
    ind = lines.index("*SCA491_RotAngle")
    nrs1_angle = float(lines[ind+1])
    ind = lines.index("*SCA491_PosX")
    nrs1_posx = float(lines[ind+1])
    ind = lines.index("*SCA491_PosY")
    nrs1_posy = float(lines[ind+1])

    # NRS2
    ind = lines.index("*SCA492_PitchX")
    nrs2_pitchx = float(lines[ind+1])
    ind = lines.index("*SCA492_PitchY")
    nrs2_pitchy = float(lines[ind+1])
    ind = lines.index("*SCA492_RotAngle")
    nrs2_angle = float(lines[ind+1])
    ind = lines.index("*SCA492_PosX")
    nrs2_posx = float(lines[ind+1])
    ind = lines.index("*SCA492_PosY")
    nrs2_posy = float(lines[ind+1])

    # NRS1 Sky to Detector
    scaling = np.array([[1/nrs1_pitchx, 0], [0, 1/nrs1_pitchy]])
    rotmat = models.Rotation2D._compute_matrix(-nrs1_angle)
    matrix = np.dot(scaling, rotmat)
    aff = models.AffineTransformation2D(matrix, name='fpa_affine_s2d')
    nrs1_sky2det = models.Shift(-nrs1_posx, name='fpa_x_s2d') & \
        models.Shift(-nrs1_posy, name='fpa_y_s2d') | aff

    # NRS1 Detector to Sky
    rotmat = models.Rotation2D._compute_matrix(-nrs1_angle)
    scaling = np.array([[nrs1_pitchx, 0], [0, nrs1_pitchy]])
    matrix = np.dot(rotmat, scaling)
    aff = models.AffineTransformation2D(matrix, name='fpa_affine_d2s')
    nrs1_det2sky = aff | models.Shift(nrs1_posx, name='fpa_x_d2s') & \
        models.Shift(nrs1_posy, name='fpa_y_d2s')

    nrs1_det2sky.inverse = nrs1_sky2det

    # NRS2 Sky to Detector
    scaling = np.array([[-1/nrs2_pitchx, 0], [0, -1/nrs2_pitchy]])
    rotmat = models.Rotation2D._compute_matrix(-nrs2_angle)
    matrix = np.dot(scaling, rotmat)
    aff = models.AffineTransformation2D(matrix, name='fpa_affine_s2d')
    nrs2_sky2det = models.Shift(-nrs2_posx, name='fpa_x_s2d') & \
        models.Shift(-nrs2_posy, name='fpa_y_s2d') | aff

    # NRS2 Detector to Sky
    rotmat = models.Rotation2D._compute_matrix(nrs2_angle)
    scaling = np.array([[-nrs2_pitchx, 0], [0, -nrs2_pitchy]])
    matrix = np.dot(scaling, rotmat)
    aff = models.AffineTransformation2D(matrix, name='fpa_affine_d2s')
    nrs2_det2sky = aff | models.Shift(nrs2_posx, name='fpa_x_d2s') &  \
        models.Shift(nrs2_posy, name='fpa_y_d2s')

    nrs2_det2sky.inverse = nrs2_sky2det

    fpa_model = FPAModel()
    fpa_model.nrs1_model = nrs1_det2sky
    fpa_model.nrs2_model = nrs2_det2sky
    fpa_model.meta.author = author
    fpa_model.meta.description = description
    fpa_model.meta.useafter = useafter
    fpa_model.meta.pedigree = "GROUND"

    return fpa_model

def create_fpa_reference(fpa_refname, out_name, author=None, description=None, useafter=None):
    with open(fpa_refname) as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]
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

    try:
        model = fpa2asdf(fpa_refname, author, description, useafter)
    except:
        raise
    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    model.history.append(entry)
    model.to_asdf(out_name)
    model.validate()


'''
if __name__ == '__main__':
    import argpars
    parser = argpars.ArgumentParser(description="Creates NIRSpec 'fpa' reference file in ASDF format.")
    parser.add_argument("fpa_file", type=str, help="FPA file.")
    parser.add_argument("output_name", type=str, help="Output file name")
    res = parser.parse_args()
    if res.output_name is None:
        output_name = "nirspec_fpa.asdf"
    else:
        output_name = res.output_name
    ref_kw = common_reference_file_keywords("FPA", "NIRSPEC FPA Description - CDP4")

    try:
        fpa2asdf(res.fpa_file, output_name, ref_kw)
    except:
        raise Exception("FPA file was not converted.")
'''
