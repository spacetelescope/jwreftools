from asdf import AsdfFile
from astropy.modeling.models import Mapping, Identity
from astropy.modeling import models
from .utils import (linear_from_pcf_det2sky, common_reference_file_keywords,
                    coeffs_from_pcf)


def slitfore2asdf():
    # fore reference file
    for i, filter in enumerate(["CLEAR", "F070LP", "F100LP", "F110W", "F140X", "F170LP", "F290LP"]):
        fore_name = "nirspec_fore_{0}.asdf".format(filter)
        ref_kw = common_reference_file_keywords(reftype="fore", title="NIRSPEC FORE Model - CDP4",
                                                description="FORE transform.",
                                                exp_type="N/A", useafter=useafter, author=author,
                                                filename=fore_name)
        ref_kw['filter'] = filter
        filename = "Fore_{0}.pcf".format(filter)
        fore_refname = os.path.join(ref_files, "CoordTransform", filename)
        
        try:
            fore2asdf(fore_refname, fore_name, ref_kw)
        except:
            raise Exception(("FORE file was not created - filter {0}".format(filter)))
            


def ifufore2asdf():
    filename = "IFU_FORE.pcf"
    ifufore_refname = os.path.join(ref_files, "CoordTransform", "IFU", filename)

    try:
        _fore2asdf(ifufore_refname, ifufore_name, ref_kw)
    except:
        print("FORE file was not created.")
        raise


def fore2asdf(pcffore, outname, ref_kw):
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

    model = (Mapping((0, 1, 2, 0, 1, 2)) | model_poly) | fore_linear

    f = AsdfFile()
    f.tree = ref_kw.copy()
    f.tree['model'] = model
    f.add_history_entry("Build 6")
    asdffile = f.write_to(outname)
    return asdffile


if __name__ == '__main__':
    import argpars
    parser = argpare.ArgumentParser(description="Creates NIRSpec 'fore' reference files in ASDF format.")
    parser.add_argument("mode", type=str, help="'ifu' or 'slit'")

    res = parser.parse_args()
    if res.mode == "ifu":
        ifufore2asdf()
    else:
        slitfore2asdf()
