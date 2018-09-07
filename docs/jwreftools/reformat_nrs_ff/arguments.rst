Arguments
=========
There is one user-callable function in the `reformat_nrs_ff` module,
`process_files`.  The script may be run by importing reformat_nrs_ff
in Python, then executing reformat_nrs_ff.process_files().  The files
to be converted have (previously, at least) been provided in several
directories, so it was intended that this script be run once for each
of these directories, using a wildcard to select all FITS files in
each directory.  The output may be written to the default directory.
For example::

    import reformat_nrs_ff
    reformat_nrs_ff.process_files("directory1/*fits", prefix="v2_")
    reformat_nrs_ff.process_files("directory2/*fits", prefix="v2_")
    reformat_nrs_ff.process_files("directory3/*fits", prefix="v2_")

process_files
-------------
This function will read and convert the specified FITS files.

*  ``pattern``

``pattern`` is a string, the full path name to one or more NIRSpec
flat-field files in FITS format.  Wildcards may be used.  These files
will be copied and reformatted; the input files themselves will not be
modified.  See also `prefix`.

* ``prefix``

For each input file name, the output file name will be constructed by
taking the basename of the input file name, and then appending the
basename to `prefix`.  Note that unless `prefix` includes a directory
name, the output file or files will be written to the default directory.
If `prefix` does include a directory name, that directory must already
exist.  The input and output names may be the same if prefix="" and the
input name includes a directory that is different from the default
directory.

* ``verbose``

If True, information will be printed to the standard output.
