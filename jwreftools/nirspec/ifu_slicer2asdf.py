import os.path
from astropy.modeling import models
from asdf import AsdfFile
from .utils import common__reference_file_keywords


def ifu_slicer2asdf(ifuslicer, outname, ref_kw):
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
    tree = ref_kw.copy()
    data = f[1].data
    header = f[1].header
    shiftx = models.Shift(header['XREF'], name='ifu_slicer_xref')
    shifty = models.Shift(header['YREF'], name='ifu_slicer_yref')
    rot = models.Rotation2D(np.rad2deg(header['ROT']), name='ifu_slicer_rot')
    model = rot | shiftx & shifty
    tree['model'] = model
    tree['data'] = f[1].data
    f.close()
    fasdf = AsdfFile()
    fasdf.tree = tree
    fasdf.add_history_entry("Build 6")
    fasdf.write_to(outname)
    return fasdf


if __name__ == '__main__':
    import argpars
    parser = argpare.ArgumentParser(description="Creates NIRSpec 'ifuslicer' reference file in ASDF format.")
    parser.add_argument("ifuslicer_file", type=str, help="IFU Slicer file")
    parser.add_argument("output_name", type=str, help="Output file name")
    res = parser.parse_args()
    if res.output_name is None:
        output_name = "nirspec_ifuslicer.asdf"
    else:
        output_name = res.output_name
    ref_kw = common_reference_file_keywords(reftype="ifuslicer", title="NIRSPEC IFU SLICER description - CDP4",
                                            description="Perfect slicer with 30 slices of size of 0.8 mm x 12 mm. No tilt.",
                                            exp_type="NRS_IFU", useafter=useafter, author=author,
                                            filename=output_name)

    try:
        ifu_slicer2asdf(res.ifuslicer_file, output_name, ref_kw)
    except:
        raise Exception("IFUSLICER file was not created.")
