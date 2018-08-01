import os.path
from . import *

model_dir = "/internal/1/astropy/jwreftools/cv3"


def generate(model_dir):
    # Create FPA file
    fpa_refname = os.path.join(model_dir, "Description", "FPA.fpa")
    create_fpa_reference(fpa_refname, "nirspec_cv3_fpa.asdf")

    # Create CAMERA file 
    camera_refname = os.path.join(model_dir, "CoordTransform", "Camera.pcf")
    create_camera_reference(camera_refname, "nirspec_cv3_camera.asdf")

    # Create COLLIMATOR file
    collimator_refname = os.path.join(model_dir, "CoordTransform", "Collimator.pcf")
    create_collimator_reference(collimator_refname, "nirspec_cv3_collimator.asdf")
    
    # Create DISPERSER files
    create_disperser_refs(model_dir)
    
    # Create FORE files
    create_fore_reference(model_dir)
    
    # Create IFUFORE file
    ifufore_refname = os.path.join(model_dir, "CoordTransform", "IFU", "IFU_FORE.pcf")
    create_ifufore_reference(ifufore_refname, "nirspec_cv3_ifufore.asdf")

    # Create OTE file
    ote_refname = os.path.join(model_dir, "CoordTransform", "OTE.pcf")
    create_ote_reference(ote_refname, "nirspec_cv3_ote.asdf")

    # Create IFUSLICER file
    ifuslicer_refname = os.path.join(model_dir, "Description", "IFU_slicer.sgd")
    create_ifuslicer_reference(ifuslicer_refname, "nirspec_cv3_ifuslicer.asdf")

    # Create IFUPOST file 
    create_ifupost_reference(model_dir, "nirspec_cv3_ifupost.asdf")

    # Create MSA file 
    msa_refname = os.path.join(model_dir, "Description", "MSA.msa")
    create_msa_reference(msa_refname, "nirspec_cv3_msa.asdf")

    # Create WAVELENGTHRANGE file
    wrange_refname = os.path.join(model_dir, "spectralconfigurations_rev1.2.txt")
    create_wavelengthrange_reference(wrange_refname, "nirspec_cv3_wavelengthrange.asdf")
