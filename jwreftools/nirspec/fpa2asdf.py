import numpy as np
from asdf import AsdfFile
from astropy.modeling import models
from .utils import common_reference_file_keywords


def fpa2asdf(fpafile, outname, ref_kw):
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
    #print('nrs2_angle', nrs2_angle, nrs2_pitchx, nrs2_pitchy, nrs2_posx, nrs2_posy)
    tree = ref_kw.copy()

    # NRS1 Sky to Detector
    scaling = np.array([[1/nrs1_pitchx, 0], [0, 1/nrs1_pitchy]])
    rotmat = models.Rotation2D._compute_matrix(-nrs1_angle)
    matrix = np.dot(scaling, rotmat)
    aff = models.AffineTransformation2D(matrix, name='fpa_affine_sky2detector')
    nrs1_sky2det = models.Shift(-nrs1_posx, name='fpa_shift_x') & \
        models.Shift(-nrs1_posy, name='fpa_shift_y') | aff

    # NRS1 Detector to Sky
    rotmat = models.Rotation2D._compute_matrix(-nrs1_angle)
    scaling = np.array([[nrs1_pitchx, 0], [0, nrs1_pitchy]])
    matrix = np.dot(rotmat, scaling)
    aff = models.AffineTransformation2D(matrix, name='fpa_affine_detector2sky')
    nrs1_det2sky = aff | models.Shift(nrs1_posx, name='fpa_shift_x_det2sky') & \
        models.Shift(nrs1_posy, name='fpa_shift_y_det2sky')

    nrs1_det2sky.inverse = nrs1_sky2det

    # NRS2 Sky to Detector
    scaling = np.array([[-1/nrs2_pitchx, 0], [0, -1/nrs2_pitchy]])
    rotmat = models.Rotation2D._compute_matrix(-nrs2_angle)
    matrix = np.dot(scaling, rotmat)
    aff = models.AffineTransformation2D(matrix, name='fpa_affine_sky2detector')
    nrs2_sky2det = models.Shift(-nrs2_posx, name='fpa_shixft_x') & \
        models.Shift(-nrs2_posy, name='fpa_shift_y') | aff
    
    # NRS2 Detector to Sky
    rotmat = models.Rotation2D._compute_matrix(nrs2_angle)
    scaling = np.array([[-nrs2_pitchx, 0], [0, -nrs2_pitchy]])
    matrix = np.dot(scaling, rotmat)
    aff = models.AffineTransformation2D(matrix, name='fpa_affine_detector2sky')
    nrs2_det2sky = aff | models.Shift(nrs2_posx, name='fpa_shift_x_det2sky') & \
        models.Shift(nrs2_posy, name='fpa_shift_y_det2sky')
    
    nrs2_det2sky.inverse = nrs2_sky2det

    tree['NRS1'] = nrs1_det2sky
    tree['NRS2'] = nrs2_det2sky
    fasdf = AsdfFile()
    fasdf.tree = tree
    fasdf.add_history_entry("Build 6")
    fasdf.write_to(outname)
    return fasdf

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

