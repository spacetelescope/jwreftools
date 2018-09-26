Description
===========
This module converts NIRSpec flat-field reference files from the IDT's
format to the format that can be submitted to CRDS for use by the
jwst/flat_field step.

Here is a summary of the differences to be expected between the input
and output files.

In the fore optics flats for MSA data, the extension names include "_Q1",
"_Q2", "_Q3", or "_Q4" (e.g. SCI_Q1), depending on the quadrant.  The
reformatted file will use EXTNAME = "SCI" (for example), together with
EXTVER = 1, 2, 3, or 4 to identify the quadrant.

The tables for "fast" variation of the flat field with wavelength have
extension names RQE, IFU, or VECTOR.  One of these tables has column name
"RQE" for the flat-field data, while all other such tables have column name
"DATA".  In the output files, the extension name for this table will be
"FAST_VARIATION", and the column name for flat-field data will be "DATA".
In some of the fast-variation tables, the wavelength and data values are 0;
in this case, wavelength and data will be changed to 1.

For fixed-slit data, separate extensions are used for the fast variation
tables, with EXTNAME set to the slit name (but not the slit names used by
DMS).  The output file will have one fast-variation table with with
multiple rows, with column names "wavelength" and "data" (the latter is
for the flat-field values), and with the DMS slit name in a column with
name "slit_name".

The DQ_DEF extensions have column names "BIT", "VAL", "NAME", and
"DESCRIPTION".  "VAL" will be renamed to "VALUE".  The string lengths for
NAME and DESCRIPTION will be assigned new values.

For SCI extensions that have a 3-D data array, each plane corresponds to
a different wavelength.  The wavelength value for each plane was written
to a keyword in the SCI extension header.  In some files the keywords
started with "FLAT", and in other files the keywords started with "PFLAT".
Following "FLAT" or "PFLAT", there was an integer index that started
with 0 in some files and 1 in other files.  In the output files, the
wavelengths will be stored in a new table extension with a single column,
"WAVELENGTH".

The flat field image data (i.e. not the tables) are in detector
orientation.  These will be changed to DMS (sky) orientation in the output.

In some files, the DQ image extension has data type float; this will be
cast to int32 for the output file.

Some required primary header keywords will be added, renamed, or modified.
INSTRUME was occasionally written INSRTUME.  USEAFTER needs to include the
time (even if it's 0), and the value was in some cases not a date at all.
EXP_TYPE will be added.  REFTYPE will be changed to FFLAT, SFLAT, or DFLAT.
If FILTER is OPAQUE, it will be changed to the value that would be expected
for a science exposure.
