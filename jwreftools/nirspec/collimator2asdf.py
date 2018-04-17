import datetime
from .utils import pcf2model#, common_reference_file_keywords)
from jwst.datamodels import CollimatorModel
from asdf.tags.core import Software, HistoryEntry

__all__ = ["create_collimator_reference", "collimator2asdf"]


def collimator2asdf(collimator_refname, author, description, useafter):
    try:
        model = pcf2model(collimator_refname, name="collimator")
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
    collimator_model.history['entries'] = [entry]

    return collimator_model


def create_collimator_reference(collimator_refname, out_name, author=None, description=None, useafter=None):
    with open(collimator_refname) as f:
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
        model = collimator2asdf(collimator_refname, author, description, useafter)
    except:
        raise
    model.to_asdf(out_name)
    new_model = CollimatorModel(out_name)
    new_model.validate()
