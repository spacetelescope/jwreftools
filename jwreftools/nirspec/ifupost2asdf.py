import datetime
import os.path, glob
from asdf.tags.core import Software, HistoryEntry
from astropy.modeling import models
from astropy.modeling.models import Mapping, Identity
from jwst.datamodels import IFUPostModel
from .utils import (common_reference_file_keywords, coeffs_from_pcf,
                    homothetic_sky2det)

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
        linear_sky2det = homothetic_sky2det(input_rot_center, rotation_angle, factors, output_rot_center)

        degree = int(lines[lines.index('*FitOrder') + 1])

        xcoeff_index = lines.index('*xForwardCoefficients 21 2')
        xlines = lines[xcoeff_index + 1: xcoeff_index + 22]
        xcoeff_forward = coeffs_from_pcf(degree, xlines)
        x_poly_forward = models.Polynomial2D(degree, name='x_ifupost_forward', **xcoeff_forward)

        ycoeff_index = lines.index('*yForwardCoefficients 21 2')
        ycoeff_forward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
        y_poly_forward = models.Polynomial2D(degree, name='y_ifupost_forward', **ycoeff_forward)

        xcoeff_index = lines.index('*xBackwardCoefficients 21 2')
        xcoeff_backward = coeffs_from_pcf(degree, lines[xcoeff_index + 1: xcoeff_index + 22])
        x_poly_backward = models.Polynomial2D(degree, name='x_ifupost_backward', **xcoeff_backward)

        ycoeff_index = lines.index('*yBackwardCoefficients 21 2')
        ycoeff_backward = coeffs_from_pcf(degree, lines[ycoeff_index + 1: ycoeff_index + 22])
        y_poly_backward = models.Polynomial2D(degree, name='y_ifupost_backward', **ycoeff_backward)

        output2poly_mapping = Identity(2, name='output_mapping_ifupost')
        output2poly_mapping.inverse = Mapping([0, 1, 0, 1])
        input2poly_mapping = Mapping([0, 1, 0, 1], name='input_mapping_ifupost')
        input2poly_mapping.inverse = Identity(2)

        model_poly = input2poly_mapping | (x_poly_forward & y_poly_forward) | output2poly_mapping

        model = linear_sky2det | model_poly
        name = "slice_{0}".format(n)
        setattr(ifupost_model, name, model)
    
    ifupost_model.meta.author = author
    ifupost_model.meta.description = description
    ifupost_model.meta.useafter = useafter
    ifupost_model.meta.pedigree = "GROUND"
    return ifupost_model



#if __name__ == '__main__':
#    import argpars
#    parser = argpare.ArgumentParser(description="Creates NIRSpec 'ifupost' reference file in ASDF format.")
#    parser.add_argument("ifupost_file_list", type=(list, str), help="IFUPOST file list or directory with IFUPOST files.")
#    parser.add_argument("output_name", type=str, help="Output file name")
#    res = parser.parse_args()
#    if res.output_name is None:
#        output_name = "nirspec_ifupost.asdf"
#    else:
#        output_name = res.output_name
#    ref_kw = common_reference_file_keywords(reftype="ifupost", title="NIRSPEC IFU-POST transforms - CDP4",
#                                           description="Cold design IFU transform, cout x+y fitted with FF/Argon/IMA exposures.",
#                                            exp_type="NRS_IFU", useafter=useafter, author=author,
#                                            filename=ifupost_name)

def create_ifupost_reference(model_dir, out_name, author=None, description=None, useafter=None):
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
    model.history = [entry]
    model.to_asdf(out_name)
