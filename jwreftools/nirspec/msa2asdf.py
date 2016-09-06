import os.path
from astropy.modeling import models
from .utils import common_reference_file_keywords
 

def msa2asdf(msafile, outname, ref_kw):
    """
    Create an asdf reference file with the MSA description.

    mas2asfdf("MSA.msa", "msa.asdf")

    Parameters
    ----------
    msafile : str
        A fits file with MSA description (MSA.msa)
    outname : str
        Name of output ASDF file.
    """
    f = fits.open(msafile)
    tree = ref_kw.copy()
    data = f[5].data # SLITS and IFU
    header = f[5].header
    shiftx = models.Shift(header['SLITXREF'], name='slit_xref')
    shifty = models.Shift(header['SLITYREF'], name='slit_yref')
    slitrot = models.Rotation2D(np.rad2deg(header['SLITROT']), name='slit_rot')

    tree[5] = {}
    tree[5]['model'] = slitrot | shiftx & shifty
    tree[5]['data'] = f[5].data
    for i in range(1, 5):
        header = f[i].header
        shiftx = models.Shift(header['QUADXREF'], name='msa_xref')
        shifty = models.Shift(header['QUADYREF'], name='msa_yref')
        slitrot = models.Rotation2D(np.rad2deg(header['QUADROT']), name='msa_rot')
        tree[i] = {}
        tree[i]['model'] = slitrot | shiftx & shifty
        tree[i]['data'] = f[i].data

    f.close()
    fasdf = AsdfFile()
    fasdf.tree = tree
    fasdf.add_history_entry("Build 6")
    fasdf.write_to(outname)
    return fasdf


if __name__ == '__main__':
    import argpars
    parser = argpars.ArgumentParser(description="Creates NIRSpec MSA reference file in ASDF format.")
    parser.add_argument("msa_file", type=str, help="MSA description file.")
    parser.add_argument("output_name", type=str, help="Output file name")
    res = parser.parse_args()
    if res.output_name is None:
        output_name = "nirspec_msa.asdf"
    else:
        output_name = res.output_name
   ref_kw = common_reference_file_keywords("MSA", "NIRSPEC MSA Description - CDP4")

   try:
       msa2asdf(res.msa_file, output_name, ref_kw)
   except:
     raise Exception("MSA file was not converted")

