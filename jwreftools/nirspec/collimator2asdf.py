import datetime
from .utils import pcf2model#, common_reference_file_keywords)
from jwst.datamodels import CollimatorModel
from asdf.tags.core import Software, HistoryEntry


def collimator2asdf(collimator_refname, author, description, useafter):
    try:
        model = pcf2model(collimator_refname)
    except:
        print("Collimator file was not converted.")
        raise
    collimator_model = CollimatorModel(model=model)
    collimator_model.meta.author = author
    collimator_model.meta.description = description
    collimator_model.meta.pedigree = "GROUND"
    collimator_model.meta.title = "NIRSPEC COLLIMATOR file"
    collimator_model.meta.useafter = useafter
    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    collimator_model.history = [entry]

    return collimator_model


#if __name__ == '__main__':
    #import argpars
    #parser = argpars.ArgumentParser(description="Creates NIRSpec 'collimator' reference file in ASDF format.")
    #parser.add_argument("collimator_file", type=str, help="Collimator file.")
    #parser.add_argument("output_name", type=str, help="Output file name")
    #res = parser.parse_args()
    #if res.output_name is None:
        #output_name = "nirspec_collimator.asdf"
    #else:
        #output_name = res.output_name

    #ref_kw = common_reference_file_keywords("COLLIMATOR", "NIRSPEC Collimator Model - CDP4")
    #collimator2asdf(collimator_file, output_name, ref_kw)

