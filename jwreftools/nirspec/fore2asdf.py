import datetime
import os.path
from asdf.tags.core import Software, HistoryEntry
from astropy.modeling.models import Mapping, Identity
from astropy.modeling import models
from .utils import linear_from_pcf_det2sky, coeffs_from_pcf

from jwst.datamodels import FOREModel, IFUFOREModel


__all__ = ["create_fore_reference", "create_ifufore_reference", "fore2asdf"]


def create_fore_reference(refdir, author=None, description=None, useafter=None):
    # fore reference file
    for i, filter in enumerate(["CLEAR", "F070LP", "F100LP", "F110W", "F140X", "F170LP", "F290LP"]):
        filename = "Fore_{0}.pcf".format(filter)
        out_name = "fore_cv3_{0}.asdf".format(filter)
        fore_refname = os.path.join(refdir, "CoordTransform", filename)
        with open(fore_refname) as f:
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
            model = fore2asdf(fore_refname)
        except:
            raise Exception(("FORE file was not created - filter {0}".format(filter)))
        fore_model = FOREModel()
        fore_model.model = model
        fore_model.meta.pedigree = 'GROUND'
        fore_model.meta.author = author
        fore_model.meta.description = description
        fore_model.meta.useafter = useafter
        fore_model.meta.instrument.filter = filter

        entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
        software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                             'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
        entry['software'] = software
        fore_model.history = [entry]
        fore_model.to_asdf(out_name)


def create_ifufore_reference(ifufore_refname, out_name, author=None, description=None, useafter=None):
    #filename = "IFU_FORE.pcf"
    with open(ifufore_refname) as f:
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
        model = fore2asdf(ifufore_refname)
    except:
        print("FORE file was not created.")
        raise
    ifufore_model = IFUFOREModel()
    ifufore_model.model = model
    ifufore_model.meta.pedigree = 'GROUND'
    ifufore_model.meta.author = author
    ifufore_model.meta.description = description
    ifufore_model.meta.useafter = useafter

    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    ifufore_model.history = [entry]
    ifufore_model.to_asdf(out_name)


def fore2asdf(pcffore):
    """
    forward direction : msa 2 ote
    backward_direction: msa 2 fpa

    """
    with open(pcffore) as f:
        lines = [l.strip() for l in f.readlines()]

    fore_det2sky = linear_from_pcf_det2sky(pcffore)
    fore_linear = fore_det2sky
    fore_linear.inverse = fore_det2sky.inverse & Identity(1)

    # compute the polynomial
    degree = int(lines[lines.index('*FitOrder') + 1])

    xcoeff_index = lines.index('*xForwardCoefficients 21 2')
    xlines = lines[xcoeff_index + 1: xcoeff_index + 22]
    xcoeff_forward = coeffs_from_pcf(degree, xlines)
    # Polynomial Correction in x
    x_poly_backward = models.Polynomial2D(degree, name="x_poly_backward", **xcoeff_forward)
    xlines_distortion = lines[xcoeff_index + 22: xcoeff_index + 43]
    xcoeff_forward_distortion = coeffs_from_pcf(degree, xlines_distortion)
    x_poly_backward_distortion = models.Polynomial2D(degree, name="x_backward_distortion",
                                                     **xcoeff_forward_distortion)
    # do chromatic correction
    # the input is Xote, Yote, lam
    model_x_backward = (Mapping((0, 1), n_inputs=3) | x_poly_backward) + \
                     ((Mapping((0,1), n_inputs=3) | x_poly_backward_distortion) * \
                     (Mapping((2,)) | Identity(1)))

    ycoeff_index = lines.index('*yForwardCoefficients 21 2')
    ycoeff_forward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
    y_poly_backward = models.Polynomial2D(degree, name="y_poly_backward",  **ycoeff_forward)

    ylines_distortion = lines[ycoeff_index + 22: ycoeff_index + 43]
    ycoeff_forward_distortion = coeffs_from_pcf(degree, ylines_distortion)
    y_poly_backward_distortion = models.Polynomial2D(degree, name="y_backward_distortion",
                                                     **ycoeff_forward_distortion)

    # do chromatic correction
    # the input is Xote, Yote, lam
    model_y_backward = (Mapping((0,1), n_inputs=3) | y_poly_backward) + \
                     ((Mapping((0, 1), n_inputs=3) | y_poly_backward_distortion) * \
                     (Mapping((2,)) | Identity(1) ))

    xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
    xcoeff_backward = coeffs_from_pcf(degree, lines[xcoeff_index + 1: xcoeff_index + 22])
    x_poly_forward = models.Polynomial2D(degree,name="x_poly_forward", **xcoeff_backward)

    xcoeff_backward_distortion = coeffs_from_pcf(degree, lines[xcoeff_index + 22: xcoeff_index + 43])
    x_poly_forward_distortion = models.Polynomial2D(degree, name="x_forward_distortion", **xcoeff_backward_distortion)

    # the chromatic correction is done here
    # the input is Xmsa, Ymsa, lam
    model_x_forward = (Mapping((0,1), n_inputs=3) | x_poly_forward) + \
                    ((Mapping((0,1), n_inputs=3) | x_poly_forward_distortion) * \
                    (Mapping((2,)) | Identity(1)))

    ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
    ycoeff_backward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
    y_poly_forward = models.Polynomial2D(degree, name="y_poly_forward",**ycoeff_backward)

    ycoeff_backward_distortion = coeffs_from_pcf(degree, lines[ycoeff_index + 22: ycoeff_index + 43])
    y_poly_forward_distortion = models.Polynomial2D(degree, name="y_forward_distortion",
                                                    **ycoeff_backward_distortion)

    # do chromatic correction
    # the input is Xmsa, Ymsa, lam
    model_y_forward = (Mapping((0,1), n_inputs=3) | y_poly_forward) + \
                    ((Mapping((0,1), n_inputs=3) | y_poly_forward_distortion) * \
                    (Mapping((2,)) | Identity(1)))

    #assign inverse transforms
    model_x = model_x_forward.copy()
    model_y = model_y_forward.copy()

    model_x.inverse = model_x_backward
    model_y.inverse = model_y_backward

    output2poly_mapping = Identity(2, name="output_mapping")
    output2poly_mapping.inverse = Mapping([0, 1, 2, 0, 1, 2])
    input2poly_mapping = Mapping([0, 1, 2, 0, 1, 2], name="input_mapping")
    input2poly_mapping.inverse = Identity(2)

    model_poly = input2poly_mapping  | (model_x & model_y) | output2poly_mapping
    model = model_poly | fore_linear

    return model
