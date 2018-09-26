import glob
import os

import numpy as np
from astropy.io import fits

__all__ = ['process_files']


# An output file name will be constructed by appending the input file name
# to a prefix; this is the default string.
PREFIX = "v5_"

# If there is no USEAFTER keyword in the primary header, a keyword will be
# created with this value.
NEW_USEAFTER = "2016-06-20T00:00:00"

# Function reformat_one_file determines the exposure type and the
# optical path part from strings in the input file name, as follows:
#  exposure type:
#    "FS"  --> "NRS_FIXEDSLIT"
#    "IFU" --> "NRS_IFU"
#    "MOS" --> "NRS_MSASPEC"
#  optical path part:
#    "fflat" --> fore optics part
#    "sflat" --> spectrograph part
#    "dflat" --> detector part

# Possible values for optical_path_part.
F_FLAT = "fflat"
S_FLAT = "sflat"
D_FLAT = "dflat"

# Mapping from the slit names used by the IDT to the slit names used by DMS.
dms_slit_names = {"SLIT_A_200_1": "S200A1",
                  "SLIT_A_200_2": "S200A2",
                  "SLIT_A_400": "S400A1",
                  "SLIT_B_200": "S200B1",
                  "SLIT_A_1600": "S1600A1"}

# For S_FLAT files for FS and IFU, the FILTER keyword is "OPAQUE".  The
# correct filter name can be determined from the FLAT[1-5] portion of the
# file name (or from the REFTYPE keyword, e.g. CAAFLAT1).
filter_names = {"FLAT1": "F100LP",
                "FLAT2": "F170LP",
                "FLAT3": "F290LP",
                "FLAT4": "F070LP",
                "FLAT5": "CLEAR"}


def process_files(pattern, prefix=PREFIX, verbose=True):
    """Convert all the FITS files in the specified directory.

    Parameters
    ----------
    pattern: str
        The path name (wildcards are permitted) to one or more NIRSpec
        flat-field files.  These files will be copied and reformatted;
        the input files themselves will not be modified.
        See also `prefix`.

    prefix: str, e.g. "v5_"
        For each input file, the output file name will be constructed by
        removing the directory (if any) and then concatenating to `prefix`.

    verbose: bool
        If True, print info.
    """

    files = glob.glob(pattern)
    if len(files) > 1:
        files.sort()
    bad_files = ["nirspec_MOS_sflat_G140H_OPAQUE_FLAT4_nrs1_f_01.01.fits"]
    for filename in files:
        if os.path.basename(filename) in bad_files:
            print("Warning:  skipping {}".format(filename), flush=True)
            continue
        reformat_one_file(filename, prefix, verbose)


def reformat_one_file(filename, prefix=PREFIX, verbose=True):
    """Convert the specified file.

    Parameters
    ----------
    filename: str
        The full name (i.e. including directory) of a NIRSpec flat-field
        reference file, in the format used by the NIRSpec IDT.  The file
        itself will not be modified; instead, a reformatted copy of the
        file will be created, using the same name but with a prefix
        prepended.  See also `prefix`.

    prefix: str, e.g. "v5_"
        The output file name will be the same as the input name, except
        that the output name will be prefixed by this string, and the output
        file will be written to the default directory (unless `prefix`
        includes a directory name), rather than to the directory containing
        the input file.

    verbose: bool
        If True, print info.
    """

    input = os.path.basename(filename)
    output = prefix + input
    if output == filename:
        raise RuntimeError("Input and output file names may not be the same.")

    if filename.find("_IFU") >= 0:
        exposure_type = "NRS_IFU"
    elif filename.find("_MOS") >= 0:
        exposure_type = "NRS_MSASPEC"
    elif filename.find("_FS") >= 0:
        exposure_type = "NRS_FIXEDSLIT"
    else:
        exposure_type = "N/A"

    filename_lc = filename.lower()
    if filename_lc.find("fflat") >= 0:
        optical_path_part = F_FLAT
    elif filename_lc.find("sflat") >= 0:
        optical_path_part = S_FLAT
    elif filename_lc.find("dflat") >= 0:
        optical_path_part = D_FLAT
    else:
        optical_path_part = None

    if verbose:
        print("filename = {}".format(filename), flush=True)
        print("exposure type = {}".format(exposure_type))
        print("optical path part = {}".format(optical_path_part))

    ifd = fits.open(filename)
    fix_primary_header(ifd, filename, NEW_USEAFTER, exposure_type,
                       optical_path_part, verbose)
    ofd = fits.HDUList(fits.PrimaryHDU(header=ifd[0].header))

    rename_extensions(ifd, exposure_type, optical_path_part, verbose)

    # sci, dq, err, etc. are lists of indices.
    (sci, dq, err, slit, fast_var, dq_def) = scan_HDU_list(ifd, verbose)

    copy_extensions_to_output(ifd, ofd, filename,
                              exposure_type, optical_path_part,
                              sci, dq, err, slit, fast_var, dq_def, verbose)

    ofd.writeto(output)
    ofd.close()
    ifd.close()


def fix_primary_header(ifd, filename, new_useafter, exposure_type,
                       optical_path_part, verbose):
    """Modify primary header keyword(s).

    If keyword INSRTUME is found in the primary header, that keyword
    will be renamed to INSTRUME.
    Keyword USEAFTER will be assigned a value.
    Keyword EXP_TYPE will be assigned a value.
    Keyword REFTYPE will be assigned a value.
    For S_FLAT files, if FILTER == "OPAQUE", a value will be assigned
    based on the "FLAT[1-5]" portion of the file name.

    Parameters
    ----------
    ifd: fits HDUList object
        The input file handle.

    filename: str
        The name of the input file.  Portions of this string may be used
        for populating the FILTER keyword.

    new_useafter: str
        This value will be used for populating USEAFTER.

    exposure_type: str
        The exposure type, e.g. "NRS_FIXEDSLIT".

    optical_path_part: str
        The component name:  "fflat", "sflat", or "dflat".

    verbose: bool
        If True, print info.
    """

    phdr = ifd[0].header
    if "INSRTUME" in phdr:
        phdr.rename_keyword("INSRTUME", "INSTRUME")
        if verbose:
            print("Keyword INSRTUME renamed to INSTRUME")

    if "USEAFTER" in phdr:
        old_useafter = phdr["useafter"]
        words = old_useafter.split("T")
        if (old_useafter.startswith("20") and len(old_useafter) == 19 and
            len(words) == 2):
            if verbose:
                print("Keyword USEAFTER was OK:  {}".format(old_useafter))
        else:
            words = old_useafter.split("-")
            if (old_useafter.startswith("20") and len(old_useafter) == 10 and
                len(words) == 3):
                    phdr["useafter"] = old_useafter + "T00:00:00"
            else:
                    phdr["useafter"] = new_useafter
            if verbose:
                print("Keyword USEAFTER changed from {} to {}"
                      .format(old_useafter, phdr["useafter"]))
    else:
        phdr["useafter"] = new_useafter
        if verbose:
            print("New keyword USEAFTER = {} was added"
                  .format(phdr["useafter"]))

    if "EXP_TYPE" in phdr:
        old_exp_type = phdr["exp_type"]
        phdr["exp_type"] = exposure_type
        if verbose:
            print("Keyword EXP_TYPE changed from {} to {}"
                  .format(old_exp_type, exposure_type))
    else:
        phdr["exp_type"] = exposure_type
        if verbose:
            print("New keyword EXP_TYPE = {} was added".format(exposure_type))

    reftype = optical_path_part.upper()
    if "REFTYPE" not in phdr:
        phdr['reftype'] = reftype
        if verbose:
            print("Keyword REFTYPE = {} was added".format(reftype))
    elif phdr['reftype'] != reftype:
        if verbose:
            print("Keyword REFTYPE changed from {} to {}"
                  .format(phdr['reftype'], reftype))
        phdr['reftype'] = reftype

    if optical_path_part == S_FLAT:
        filter = phdr['filter']
        if filter == "OPAQUE":
            words = filename.split("_")
            foundit = False
            for word in words:
                if word.startswith("FLAT"):
                    foundit = True
                    flat_id = word
                    break
            if foundit:
                new_filter = filter_names[flat_id]
                phdr['filter'] = new_filter
                print("Keyword FILTER changed from {} to {}"
                      .format(filter, new_filter))
            else:
                print("Warning:  Can't populate FILTER keyword.")


def rename_extensions(ifd, exposure_type, optical_path_part, verbose):
    """Rename some extensions.

    Parameters
    ----------
    ifd: fits HDUList object
        The input file handle.

    exposure_type: str
        The exposure type, e.g. "NRS_FIXEDSLIT".

    optical_path_part: str
        The component name:  "fflat", "sflat", or "dflat".

    verbose: bool
        If True, print info.
    """

    n_hdu = len(ifd)

    if exposure_type == "NRS_MSASPEC" and optical_path_part == F_FLAT:
        # Change extname = SCI_Q<k> to extname = SCI, extver = <k>, etc.
        # For each extension (i.e. not including the primary HDU) ...
        for i in range(1, n_hdu):
            ifd[i].header["extver"] = 1                 # default
            if ifd[i].header["extname"].upper() == "DQ_DEF":
                continue
            for k in range(1, 5):
                suffix = "Q{}".format(k)
                extname = ifd[i].header["extname"].upper()
                if extname.endswith(suffix):
                    if extname.startswith("SCI"):
                        ifd[i].header["extname"] = "SCI"
                        ifd[i].header["extver"] = k
                    elif extname.startswith("ERR"):
                        ifd[i].header["extname"] = "ERR"
                        ifd[i].header["extver"] = k
                    elif extname.startswith("DQ"):
                        ifd[i].header["extname"] = "DQ"
                        ifd[i].header["extver"] = k
                    elif extname.startswith("Q"):
                        ifd[i].header["extname"] = "FAST_VARIATION"
                        ifd[i].header["extver"] = k
                    if verbose:
                        print("{}:  {} --> {}, {}".format(i, extname,
                            ifd[i].header["extname"], ifd[i].header["extver"]))

    for i in range(1, n_hdu):
        extname = ifd[i].header["extname"].upper()
        if extname == "RQE" or extname == "IFU" or extname == "VECTOR":
            ifd[i].header["extname"] = "FAST_VARIATION"
            ifd[i].header["extver"] = 1
            if verbose:
                print("{} {} --> {} {}".format(i, extname,
                    ifd[i].header["extname"], ifd[i].header["extver"]))
        if ifd[i].header.get("extver", -999) == -999:
            ifd[i].header["extver"] = 1


def scan_HDU_list(ifd, verbose):
    """Find various extensions in the input file.

    Parameters
    ----------
    ifd: fits HDUList object
        The input file handle.

    verbose: bool
        If True, print info.

    Returns
    -------
    (sci, dq, err, slit, fast_var, dq_def): tuple of lists
        sci: list
            Indices of SCI extension(s) in the input file.
        dq: list
            Indices of DQ extension(s) in the input file.
        err: list
            Indices of ERR extension(s) in the input file.
        slit: list
            Indices of fast-variation extensions for fixed slits in the input
            file.  It will not be the case that both `slit` and `fast_var`
            are populated; both may be empty.
        fast_var: list
            Indices of fast-variation extension(s) in the input file.
        dq_def: list
            Index of the DQ_DEF extension in the input file.
    """

    n_hdu = len(ifd)

    sci = []
    dq = []
    err = []
    slit = []
    fast_var = []
    dq_def = []
    for i in range(1, n_hdu):
        extname = ifd[i].header.get("extname", "").lower()
        if extname.startswith("sci"):
            sci.append(i)
        elif extname == ("dq_def"):
            dq_def.append(i)
        elif extname.startswith("dq"):
            dq.append(i)
        elif extname.startswith("err"):
            err.append(i)
        elif extname.startswith("slit"):
            slit.append(i)
        elif extname == ("fast_variation"):
            fast_var.append(i)
        else:
            print("Note:  extension {} will be ignored, extname = {}"
                  .format(i, extname))

    if verbose:
        print("sci = {}".format(sci))
        print("dq = {}".format(dq))
        print("err = {}".format(err))
        print("slit = {}".format(slit))
        print("fast_var = {}".format(fast_var))
        print("dq_def = {}".format(dq_def), flush=True)

    return (sci, dq, err, slit, fast_var, dq_def)


def copy_extensions_to_output(ifd, ofd, filename,
                              exposure_type, optical_path_part,
                              sci, dq, err, slit, fast_var, dq_def, verbose):
    """Copy all the extensions that we want to keep.

    Parameters
    ----------
    ifd: fits HDUList object
        The input file handle.

    ofd: fits HDUList object
        The output file handle.

    filename: str
        The name of the input file.  This will be used for determining
        which detector was used, if there is no DETECTOR keyword.

    exposure_type: str
        The exposure type, e.g. "NRS_FIXEDSLIT".

    optical_path_part: str
        The component name:  "fflat", "sflat", or "dflat".

    sci: list
        Indices of SCI extension(s) in the input file.

    dq: list
        Indices of DQ extension(s) in the input file.

    err: list
        Indices of ERR extension(s) in the input file.

    slit: list
        Indices of fast-variation extensions for fixed slits in the input
        file.  It should not be the case that both `slit` and `fast_var`
        are populated; both may be empty.

    fast_var: list
        Indices of fast-variation extension(s) in the input file.

    dq_def: list
        Index of the DQ_DEF extension in the input file.

    verbose: bool
        If True, print info.
    """

    try:
        detector = ifd[0].header['detector']
    except KeyError:
        print("Warning:  DETECTOR keyword not found; using the file name.")
        filename_uc = filename.upper()
        if filename_uc.find("_NRS1") >= 0:
            detector = "NRS1"
        elif filename_uc.find("_NRS2") >= 0:
            detector = "NRS2"
        else:
            detector = "N/A"
            if verbose and optical_path_part != F_FLAT:
                print("Warning:  Could not determine which detector was used.")
    if verbose:
        print("detector = {}".format(detector), flush=True)

    if exposure_type == "NRS_MSASPEC" and optical_path_part == F_FLAT:
        shape = ifd[sci[0]].data.shape
        # The DQ arrays in the input are all ones; replace with zeros.
        dq_data = np.zeros(shape, dtype=np.uint32)
        for k in range(4):
            # Index of (SCI,k) extension.
            i = sci[k]
            imhdu = fits.ImageHDU(data=ifd[i].data.astype(np.float32),
                                  header=ifd[i].header)
            ofd.append(imhdu)
            if len(ifd[i].data.shape) > 2:
                wl_hdu = create_wavelength_table(ifd[i].header,
                                                 ifd[i].data.shape[0], k + 1)
            else:
                wl_hdu = None
            # Index of (DQ,k) extension.
            i = dq[k]
            imhdu = fits.ImageHDU(data=dq_data, header=ifd[i].header)
            ofd.append(imhdu)
            # Index of (ERR,k) extension.
            i = err[k]
            imhdu = fits.ImageHDU(data=ifd[i].data.astype(np.float32),
                                  header=ifd[i].header)
            ofd.append(imhdu)
            if wl_hdu is not None:
                ofd.append(wl_hdu)
                wl_hdu = None
            # Index of the fast-variation table extension.
            i = fast_var[k]
            fv_hdu = reformat_fast_var_table(ifd, [i], k + 1)
            ofd.append(fv_hdu)
    else:
        wl_hdu = None                           # initial value
        # Index of SCI extension.
        if len(sci) > 0:
            i = sci[0]
            sci_data = ifd[i].data
            ifd[i].data[...] = flip_sci(sci_data, detector).copy()
            ofd.append(ifd[i])
            del sci_data
            if len(ifd[i].data.shape) > 2:
                wl_hdu = create_wavelength_table(ifd[i].header,
                                                 ifd[i].data.shape[0], 1)
        # Index of DQ extension.
        if len(dq) > 0:
            i = dq[0]
            dq_data = ifd[i].data.astype(np.uint32)
            mod_data = flip_sci(dq_data, detector)
            imhdu = fits.ImageHDU(data=mod_data, header=ifd[i].header)
            del dq_data, mod_data
            ofd.append(imhdu)
        # Index of ERR extension.
        if len(err) > 0:
            i = err[0]
            err_data = ifd[i].data
            ifd[i].data[...] = flip_sci(err_data, detector).copy()
            ofd.append(ifd[i])
            del err_data
        if wl_hdu is not None:
            ofd.append(wl_hdu)
            del wl_hdu
        if len(slit) > 0:
            fv_hdu = reformat_fast_var_table(ifd, slit, 1)
            ofd.append(fv_hdu)
        if len(fast_var) > 0:
            fv_hdu = reformat_fast_var_table(ifd, fast_var, 1)
            ofd.append(fv_hdu)

    if len(dq_def) > 0:
        i = dq_def[0]
        new_dq_def = rename_dq_def_cols(ifd[i], verbose)
        ofd.append(new_dq_def)


def flip_sci(sci_data, detector):
    """Transpose and (if NRS2) rotate 180 degrees.

    Parameters
    ----------
    sci_data: numpy array
        An input data array, e.g. SCI or DQ.

    detector: str
        The detector name, either "NRS1" or "NRS2".

    Returns
    -------
    mod_data: numpy array
        The same values as in `sci_data`, but transposed and possibly
        rotated 180 degrees.
    """

    shape = sci_data.shape
    naxis = len(shape)
    if naxis < 2 or naxis > 3:
        print("Warning:  SCI data have dimension {}".format(naxis), flush=True)
        return sci_data

    if detector.upper() == "NRS1":
        if naxis == 2:
            mod_data = np.swapaxes(sci_data, 0, 1)
        elif naxis == 3:
            mod_data = np.swapaxes(sci_data, 1, 2)
    elif detector.upper() == "NRS2":
        if naxis == 2:
            mod_data = np.swapaxes(sci_data, 0, 1)[::-1, ::-1]
        elif naxis == 3:
            mod_data = np.swapaxes(sci_data, 1, 2)[:, ::-1, ::-1]

    return mod_data


def rename_dq_def_cols(input_dq_def, verbose):
    """Rename column name VAL to VALUE; change string lengths.

    Parameters
    ----------
    input_dq_def: fits HDU object
        A dq_def table

    verbose: bool
        If True, print old and new column names.

    Returns
    -------
    new_dq_def
        A copy of input_dq_def, but with column VAL renamed to Value.
    """

    for name in input_dq_def.data.names:
        if name.lower() == "value":
            if verbose:
                print("The input DQ_DEF table is OK as it is.")
            return input_dq_def

    if verbose:
        print("Input DQ_DEF: {}".format(input_dq_def.data.names))
    nrows = len(input_dq_def.data)
    col = []
    col.append(fits.Column(name="Bit", format="J"))
    col.append(fits.Column(name="Value", format="J"))
    col.append(fits.Column(name="Name", format="30A"))
    col.append(fits.Column(name="Description", format="60A"))
    # xxx col.append(fits.Column(name="BIT", format="J"))
    # xxx col.append(fits.Column(name="VALUE", format="J"))
    # xxx col.append(fits.Column(name="NAME", format="30A"))
    # xxx col.append(fits.Column(name="DESCRIPTION", format="60A"))
    cd = fits.ColDefs(col)
    new_dq_def = fits.BinTableHDU.from_columns(cd, nrows=nrows)
    new_dq_def.header["extname"] = "DQ_DEF"
    new_dq_def.header["extver"] = 1
    new_dq_def.data["BIT"][:] = input_dq_def.data["BIT"].copy()
    new_dq_def.data["VALUE"][:] = input_dq_def.data["VAL"].copy()
    new_dq_def.data["NAME"][:] = input_dq_def.data["NAME"].copy()
    new_dq_def.data["DESCRIPTION"][:] = input_dq_def.data["DESCRIPTION"].copy()
    if verbose:
        print("New DQ_DEF:   {}".format(new_dq_def.data.names), flush=True)

    return new_dq_def


def create_wavelength_table(hdr, nrows_data, extver=1):
    """Create a bintable HDU of wavelengths from header keywords.

    Parameters
    ----------
    hdr: fits Header object
        Header for a SCI extension.

    nrows_data: int
        Length of the first axis of the SCI data array.  This is the
        number of wavelength values that we expect to read from header
        keywords, and the number of rows that we will use when creating
        the table of wavelengths.  It is an error if there are not at
        least this many header keywords giving values of wavelength.

    Returns
    -------
    wl_hdu: fits BinTableHDU object
        A new HDU containing a table of wavelengths, copied from `hdr`.
    """

    k = hdr["flat*"]
    if len(k) < 1:
        k = hdr["pflat*"]
    if len(k) < 1:
        raise RuntimeError("Wavelength keywords not found.")

    nrows = nrows_data
    if len(k) != nrows_data:
        print("Note:  Data shape implies nrows = {}, "
              "but there are {} keywords for wavelength:"
              .format(nrows_data, len(k)))
        for key in k:
            print("    {}".format(key))
        if len(k) < nrows_data:
            raise RuntimeError("There must be a wavelength keyword for every "
                               "plane of the SCI data array.")
        nrows = min(nrows_data, len(k))

    wavelength = np.zeros(nrows, dtype=np.float32)
    last_wl = -1.
    in_order = True                             # initial value
    for (i, key) in enumerate(k):
        wavelength[i] = hdr[key]
        if wavelength[i] <= last_wl:
            in_order = False
        last_wl = wavelength[i]
    if not in_order:
        print("Note:  wavelengths (from keywords) are not "
              "in increasing order.", flush=True)

    col = [fits.Column(name="wavelength", format="E", unit="micrometer",
                       array=wavelength)]
    cd = fits.ColDefs(col)
    wl_hdu = fits.BinTableHDU.from_columns(cd)
    wl_hdu.header["extname"] = "WAVELENGTH"
    wl_hdu.header["extver"] = extver

    return wl_hdu


def reformat_fast_var_table(ifd, i_list, extver):
    """

    Parameters
    ----------
    ifd: fits HDUList object
        The input file handle.

    i_list: list of int
        The elements are the extension numbers of the (one or more)
        fast-variation tables.  If there are more than one, they will be
        consolidated into one output table with one row for each input
        fast-variation table.

    extver: int
        For MSA data, this identifies which of the four micro-shutter
        arrays the flat-field data are used for.  For other data, this
        will have the value 1.

    Returns
    -------
    fv_hdu: fits BinTableHDU object
        A new HDU containing a table of wavelengths and corresponding
        flat-field values.  There will be a column with name "slit_name"
        that can have value "ANY", or for fixed-slit data it will be the
        slit name (converted from the EXTNAME of the input file to the
        DMS slit name).
    """

    nrows = len(i_list)

    max_nelem = 0
    for i in i_list:
        input_wl = ifd[i].data.field("wavelength")
        max_nelem = max(max_nelem, len(input_wl))
    if max_nelem == 0:
        max_nelem = 1

    if len(i_list) > 0:
        i = i_list[0]
        fixed_slit = ifd[i].header.get("extname", "missing").startswith("SLIT")
    else:
        fixed_slit = False

    col = []
    col.append(fits.Column(name="slit_name", format="15A"))
    col.append(fits.Column(name="nelem", format="1J"))
    col.append(fits.Column(name="wavelength", format="{}E".format(max_nelem),
                           unit="micrometer"))
    col.append(fits.Column(name="data", format="{}E".format(max_nelem)))
    cd = fits.ColDefs(col)
    fv_hdu = fits.BinTableHDU.from_columns(cd, nrows=nrows)
    fv_hdu.header["extname"] = "FAST_VARIATION"
    fv_hdu.header["extver"] = extver

    row = 0
    for i in i_list:
        input_wl = ifd[i].data.field("wavelength")
        try:
            input_data = ifd[i].data.field("data")
        except KeyError:
            input_data = ifd[i].data.field("RQE")
        if np.nanmax(input_wl) <= 0.:
            print("Wavelength was {}, changed to 1.".format(input_wl[0]),
                  flush=True)
            input_wl[:] = 1.
            if max_nelem > 1:
                print("Warning:  but max_nelem = {}".format(max_nelem),
                      flush=True)
        if np.nanmax(input_data) <= 0.:
            print("Flat was {}, changed to 1.".format(input_data[0]),
                  flush=True)
            input_data[:] = 1.
            if max_nelem > 1:
                print("Warning:  but max_nelem = {}".format(max_nelem),
                      flush=True)

        nelem = len(input_wl)
        if fixed_slit:
            fv_hdu.data.field("slit_name")[row] = \
                dms_slit_names[ifd[i].header["extname"].upper()]
        else:
            fv_hdu.data.field("slit_name")[row] = "ANY"

        fv_hdu.data.field("nelem")[row] = nelem
        if max_nelem == 1:
            fv_hdu.data.field("wavelength")[:] = input_wl.copy()
            fv_hdu.data.field("data")[:] = input_data.copy()
        elif nelem > 1:
            fv_hdu.data.field("wavelength")[row, :nelem] = input_wl.copy()
            fv_hdu.data.field("wavelength")[row, nelem:] = 999.
            fv_hdu.data.field("data")[row, :nelem] = input_data.copy()
            fv_hdu.data.field("data")[row, nelem:] = 1.
        row += 1

    return fv_hdu
