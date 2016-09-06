from .utils import pcf2asdf, common_reference_file_keywords


def collimator2asdf(collimator_refname, output_name, ref_kw):
    try:
        pcf2asdf(collimator_refname, output_name_name, ref_kw)
    except:
        raise Exception("Collimator file was not converted.")


if __name__ == '__main__':
    import argpars
    parser = argpars.ArgumentParser(description="Creates NIRSpec 'collimator' reference file in ASDF format.")
    parser.add_argument("collimator_file", type=str, help="Collimator file.")
    parser.add_argument("output_name", type=str, help="Output file name")
    res = parser.parse_args()
    if res.output_name is None:
        output_name = "nirspec_collimator.asdf"
    else:
        output_name = res.output_name

    ref_kw = common_reference_file_keywords("COLLIMATOR", "NIRSPEC Collimator Model - CDP4")
    collimator2asdf(collimator_file, output_name, ref_kw)

