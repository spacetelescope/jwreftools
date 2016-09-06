from asdf import AsdfFile
from astropy.modeling import models
from astropy.modeling.momdels import Mapping, Identity
from .utils import homothetic_det2sky, common_reference_file_keywords


def ote2asdf(otepcf, outname, ref_kw):
    """
    ref_kw = common_reference_file_keywords('OTE', 'NIRSPEC OTE transform - CDP4')

    ote2asdf('Model/Ref_Files/CoordTransform/OTE.pcf', 'jwst_nirspec_ote_0001.asdf', ref_kw)
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


    f = AsdfFile()
    f.tree = ref_kw.copy()
    f.tree['model'] = model
    f.add_history_entry("Build 6")
    f.write_to(outname)
    return model_poly, mlinear

if __name__ == '__main__':
    import argpars
    parser = argpars.ArgumentParser(description="Creates NIRSpec 'ote' reference file in ASDF format.")
    parser.add_argument("ote_file", type=str, help="OTE file.")
    parser.add_argument("output_name", type=str, help="Output file name")
    res = parser.parse_args()
    if res.output_name is None:
        output_name = "nirspec_ote.asdf"
    else:
        output_name = res.output_name
    ref_kw = common_reference_file_keywords("OTE", "NIRSPEC OTE Model - CDP4")

    try:
        ote2asdf(ote_refname, ote_name, ref_kw)
    except:
        raise Exception("OTE file was not converted.")
