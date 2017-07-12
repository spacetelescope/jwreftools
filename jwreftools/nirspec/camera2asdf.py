import datetime
from .utils import pcf2model
from jwst.datamodels import CameraModel
from asdf.tags.core import Software, HistoryEntry

__all__ = ["create_camera_reference", "camera2asdf"]

def camera2asdf(camera_refname, author, description, useafter):
    try:
        model = pcf2model(camera_refname, name='camera')
    except:
        print("Camera file was not converted.")
        raise
    camera_model = CameraModel(model=model)
    camera_model.meta.author = author
    camera_model.meta.description = description
    camera_model.meta.pedigree = "GROUND"
    camera_model.meta.title = "NIRSPEC CAMERA file"
    camera_model.meta.useafter = useafter
    entry = HistoryEntry({'description': "New version created from CV3 with updated file structure", 'time': datetime.datetime.utcnow()})
    software = Software({'name': 'jwstreftools', 'author': 'N.Dencheva',
                         'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7.1"})
    entry['software'] = software
    camera_model.history = [entry]

    return camera_model

def create_camera_reference(camera_refname, out_name, author=None, description=None, useafter=None):
    with open(camera_refname) as f:
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
        model = camera2asdf(camera_refname, author, description, useafter)
    except:
        raise
    model.to_asdf(out_name)
    new_model = CameraModel(out_name)
    new_model.validate()

#if __name__ == '__main__':
    #import argpars
    #parser = argpars.ArgumentParser(description="Creates NIRSpec 'camera' reference file in ASDF format.")
    #parser.add_argument("camera_file", type=str, help="Camera file.")
    #parser.add_argument("output_name", type=str, help="Output file name")
    #res = parser.parse_args()
    #if res.output_name is None:
        #output_name = "nirspec_camera.asdf"
    #else:
        #output_name = res.output_name

    #ref_kw = common_reference_file_keywords("CAMERA", "NIRSPEC Camera Model - CDP4")

    #camera2asdf(camera_file, output_name, ref_kw)
