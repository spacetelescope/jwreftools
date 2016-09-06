from .utils import (pcf2asdf, common_reference_file_keywords)


def camera2asdf(camera_refname, output_name, ref_kw):
    try:
        pcf2asdf(camera_refname, output_name, ref_kw)
    except:
        print("Camera file was not converted.")
        raise


if __name__ == '__main__':
    import argpars
    parser = argpars.ArgumentParser(description="Creates NIRSpec 'camera' reference file in ASDF format.")
    parser.add_argument("camera_file", type=str, help="Camera file.")
    parser.add_argument("output_name", type=str, help="Output file name")
    res = parser.parse_args()
    if res.output_name is None:
        output_name = "nirspec_camera.asdf"
    else:
        output_name = res.output_name

    ref_kw = common_reference_file_keywords("CAMERA", "NIRSPEC Camera Model - CDP4")

    camera2asdf(camera_file, output_name, ref_kw)
