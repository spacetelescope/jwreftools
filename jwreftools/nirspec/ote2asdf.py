import datetime
from jwst.datamodels import OTEModel
from asdf.tags.core import Software, HistoryEntry
from astropy.modeling import models
from astropy.modeling.models import Mapping, Identity
from .utils import homothetic_det2sky, coeffs_from_pcf

__all__ = ["create_ote_reference", "ote2asdf"]


def ote2asdf(otepcf, author, description, useafter):
    """
    """
    with open(otepcf) as f:
        lines = [l.strip() for l in f.readlines()]

    factors = lines[lines.index('*Factor 2 1') + 1].split()
    # this corresponds to modeling Rotation direction as is
    rotation_angle = float(lines[lines.index('*Rotation') + 1])
    input_rot_center = lines[lines.index('*InputRotationCentre 2 1') + 1].split()
    output_rot_center = lines[lines.index('*OutputRotationCentre 2 1') + 1].split()

    mlinear = homothetic_det2sky(input_rot_center, rotation_angle, factors, output_rot_center, name="ote")

    degree = int(lines[lines.index('*FitOrder') + 1])

    xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
    xlines = lines[xcoeff_index + 1].split('\t')
    xcoeff_backward = coeffs_from_pcf(degree, xlines)
    x_poly_forward = models.Polynomial2D(degree, name='x_ote_forward', **xcoeff_backward)

    xcoeff_index = lines.index('*xForwardCoefficients 21 2')
    xlines = lines[xcoeff_index + 1].split('\t')
    xcoeff_forward = coeffs_from_pcf(degree, xlines)
    x_poly_backward = models.Polynomial2D(degree, name='x_ote_backward', **xcoeff_forward)

    ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
    ylines = lines[ycoeff_index + 1].split('\t')
    ycoeff_backward = coeffs_from_pcf(degree, ylines)
    y_poly_forward = models.Polynomial2D(degree, name='y_ote_forward', **ycoeff_backward)

    ycoeff_index = lines.index('*yForwardCoefficients 21 2')
    ylines = lines[ycoeff_index + 1].split('\t')
    ycoeff_forward = coeffs_from_pcf(degree, ylines)
    y_poly_backward = models.Polynomial2D(degree, name='y_ote_backward', **ycoeff_forward)

    x_poly_forward.inverse = x_poly_backward
    y_poly_forward.inverse = y_poly_backward

    output2poly_mapping = Identity(2, name='output_mapping_ote')
    output2poly_mapping.inverse = Mapping([0, 1, 0, 1])
    input2poly_mapping = Mapping([0, 1, 0, 1], name='input_mapping_ote')
    input2poly_mapping.inverse = Identity(2)

    model_poly = input2poly_mapping | (x_poly_forward & y_poly_forward) | output2poly_mapping

    model = model_poly | mlinear
    ote_model = OTEModel()
    ote_model.model = model
    ote_model.meta.author = author
    ote_model.meta.description = description
    ote_model.meta.useafter = useafter
    ote_model.meta.pedigree = "GROUND"
    return ote_model


def create_ote_reference(ote_file, output_name, author=None, description=None,
                         useafter="2016-03-01T09:08:05"):
    f = open(ote_file)
    lines = [l.strip() for l in f.readlines()]
    f.close()
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
        model = ote2asdf(ote_file, author, description, useafter)
    except:
        raise Exception("OTE file was not converted.")
    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    model.history = [entry]
    model.to_asdf(output_name)
