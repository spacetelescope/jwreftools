import datetime
import os.path, glob
from asdf.tags.core import Software, HistoryEntry
from astropy.modeling import models
from astropy.modeling.models import Mapping, Identity
from jwst.datamodels import IFUPostModel
from .utils import coeffs_from_pcf, homothetic_sky2det

__all__ = ["create_ifupost_reference", "ifupost2asdf"]


def ifupost2asdf(ifupost_files, author, description, useafter):
    """
    Create a reference file of type ``ifupost`` .

    Combines all IDT ``IFU-POST`` reference files in one ASDF file.

    forward direction : MSA to Collimator
    backward_direction: Collimator to MSA

    Parameters
    ----------
    ifupost_files : list
        Names of all ``IFU-POST`` IDT reference files
    outname : str
        Name of output ``ASDF`` file
    """
    ifupost_model = IFUPostModel()
    for fifu in ifupost_files:
        fifu_name = os.path.split(fifu)[1]
        n = int((fifu_name.split('IFU-POST_')[1]).split('.pcf')[0])

        with open(fifu) as f:
            lines = [l.strip() for l in f.readlines()]
        factors = lines[lines.index('*Factor 2') + 1].split()
        rotation_angle = float(lines[lines.index('*Rotation') + 1])
        input_rot_center = lines[lines.index('*InputRotationCentre 2') + 1].split()
        output_rot_center = lines[lines.index('*OutputRotationCentre 2') + 1].split()
        linear_sky2det = homothetic_sky2det(input_rot_center, rotation_angle, factors,
                                            output_rot_center, name='ifupost')

        degree = int(lines[lines.index('*FitOrder') + 1])

        xcoeff_index = lines.index('*xForwardCoefficients 21 2')
        xlines = lines[xcoeff_index + 1: xcoeff_index + 22]
        xcoeff_forward = coeffs_from_pcf(degree, xlines)
        x_poly_forward = models.Polynomial2D(degree, name='ifupost_x_forw', **xcoeff_forward)
        xlines_distortion = lines[xcoeff_index + 22: xcoeff_index + 43]
        xcoeff_forward_distortion = coeffs_from_pcf(degree, xlines_distortion)
        x_poly_forward_distortion = models.Polynomial2D(degree, name="ifupost_x_forwdist",
                                                         **xcoeff_forward_distortion)

        ycoeff_index = lines.index('*yForwardCoefficients 21 2')
        ycoeff_forward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
        y_poly_forward = models.Polynomial2D(degree, name='ifupost_y_forw', **ycoeff_forward)
        ylines_distortion = lines[ycoeff_index + 22: ycoeff_index + 43]
        ycoeff_forward_distortion = coeffs_from_pcf(degree, ylines_distortion)
        y_poly_forward_distortion = models.Polynomial2D(degree, name="ifupost_y_forwdist",
                                                     **ycoeff_forward_distortion)

        xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
        xcoeff_backward = coeffs_from_pcf(degree, lines[xcoeff_index + 1: xcoeff_index + 22])
        x_poly_backward = models.Polynomial2D(degree, name='ifupost_x_back', **xcoeff_backward)

        ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
        ycoeff_backward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
        y_poly_backward = models.Polynomial2D(degree, name='ifupost_y_back', **ycoeff_backward)

        x_poly_forward.inverse = x_poly_backward
        y_poly_forward.inverse = y_poly_backward

        output2poly_mapping = Identity(2, name='ifupost_outmap')
        output2poly_mapping.inverse = Mapping([0, 1, 2, 0, 1, 2])
        input2poly_mapping = Mapping([0, 1, 2, 0, 1, 2], name='ifupost_inmap')
        input2poly_mapping.inverse = Identity(2)

        model = {'linear': linear_sky2det,
                 'xpoly': x_poly_forward,
                 'xpoly_distortion': x_poly_forward_distortion,
                 'ypoly': y_poly_forward,
                 'ypoly_distortion': y_poly_forward_distortion
                 }
        name = "slice_{0}".format(n)
        setattr(ifupost_model, name, model)

    ifupost_model.meta.author = author
    ifupost_model.meta.description = description
    ifupost_model.meta.useafter = useafter
    ifupost_model.meta.pedigree = "GROUND"
    return ifupost_model


def create_ifupost_reference(model_dir, out_name, author=None, description=None, useafter=None):
    """
    Create the IFUPOST reference.

    Parameters
    ----------
    model_dir : str
        Directory with the model. it should contain a
        subdirectory ``CoordTransform``.
    out_name : str
        Name for the reference file.
    author : str
        Author field.
    description : str
        Consice description of the file.
    useafter : str
        A useafter date in ISO format.
    """
    model_dir = os.path.join(model_dir, "CoordTransform", "IFU")
    ifupost_list = glob.glob(model_dir + '/IFU-POST*')
    f = open(ifupost_list[0])
    lines = f.readlines()
    f.close()
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
        model = ifupost2asdf(ifupost_list, author, description, useafter)
    except:
        raise Exception("IFUPOST file was not created.")
    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    model.history.append(entry)
    model.to_asdf(out_name)
    model.validate()
