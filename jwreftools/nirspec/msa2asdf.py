import datetime
import os.path
import numpy as np
from asdf.tags.core import Software, HistoryEntry
from astropy.io import fits
from astropy.modeling import models
from jwst.datamodels import MSAModel

__all__ = ["create_msa_reference", "msa2asdf"]

def msa2asdf(msafile, author, description, useafter):
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
    data = f[5].data # SLITS and IFU
    header = f[5].header
    shiftx = models.Shift(header['SLITXREF'], name='msa_slit_x')
    shifty = models.Shift(header['SLITYREF'], name='msa_slit_y')
    slitrot = models.Rotation2D(np.rad2deg(header['SLITROT']), name='msa_slit_rot')
    msa_model = MSAModel()
    msa_model.Q5.model = slitrot | shiftx & shifty
    msa_model.Q5.data = f[5].data
    for i in range(1, 5):
        header = f[i].header
        shiftx = models.Shift(header['QUADXREF'], name='msa_Q{0}_x'.format(i))
        shifty = models.Shift(header['QUADYREF'], name='msa_Q{0}_y'.format(i))
        slitrot = models.Rotation2D(np.rad2deg(header['QUADROT']), name='msa_Q{0}_rot'.format(i))
        model = slitrot | shiftx & shifty
        data = f[i].data
        name = "Q{0}".format(i)
        setattr(msa_model, name, {'model': model, 'data': data})
    f.close()
    msa_model.meta.author = author
    msa_model.meta.description = description
    msa_model.meta.useafter = useafter
    msa_model.meta.pedigree = 'GROUND'
    return msa_model


#if __name__ == '__main__':
#    import argpars
#    parser = argpars.ArgumentParser(description="Creates NIRSpec MSA reference file in ASDF format.")
#    parser.add_argument("msa_file", type=str, help="MSA description file.")
#    parser.add_argument("output_name", type=str, help="Output file name")
#    res = parser.parse_args()
#    if res.output_name is None:
#        output_name = "nirspec_msa.asdf"
#    else:
#        output_name = res.output_name
#    ref_kw = common_reference_file_keywords("MSA", "NIRSPEC MSA Description - CDP4")
def create_msa_reference(msa_file, output_name, author=None, description=None, useafter=None):
    f = fits.open(msa_file)
    auth = f[0].header['AUTHOR']
    descrip = f[0].header['DESCR']
    date = f[0].header['DATE']
    f.close()
    if author is None:
        author = auth
    if description is None:
        description = descrip
    if useafter is None:
        useafter = date
    try:
        msa_model = msa2asdf(msa_file, author, description, useafter)
    except:
        raise Exception("MSA file was not converted")
    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    msa_model.history.append(entry)
    msa_model.to_asdf(output_name)
    msa_model.validate()

