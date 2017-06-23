import os.path
import datetime
from astropy.modeling import models

from .utils import pcf2model#, common_reference_file_keywords)
from jwst.datamodels import DisperserModel
from asdf.tags.core import Software, HistoryEntry



def disperser2asdf(disfile, tiltyfile, tiltxfile, author, description, useafter):
    """
    Create a NIRSPEC disperser reference file in ASDF format.

    Combine information stored in disperser_G?.dis and disperser_G?_TiltY.gtp
    files delievred by the IDT.

    disperser2asdf("disperser_G140H.dis", "disperser_G140H_TiltY.gtp", "disperserG140H.asdf")

    Parameters
    ----------
    disfile : list or str
        A list of .dis files or a wild card string (*.dis).
    tiltyfile : str
        File with tilt_Y data, e.g. disperser_G395H_TiltY.gtp.
    outname : str
        Name of output ASDF file.

    Returns
    -------
    fasdf : asdf.AsdfFile

    """

    disperser = disfile.split('.dis')[0].split('_')[1]
    with open(disfile) as f:
        lines=[l.strip() for l in f.readlines()]

    try:
        ind = lines.index('*TYPE')
        disperser_type = (lines[ind + 1]).lower()
    except ValueError:
        raise ValueError("Unknown disperser type in {0}".format(disfile))

    if disperser_type == 'gratingdata':
        d = dict.fromkeys(['groove_density', 'theta_z', 'theta_y', 'theta_x', 'tilt_y'])
    elif disperser_type == 'prismdata':
        d = dict.fromkeys(['tref', 'pref', 'angle', 'lcoef', 'kcoef', 'tcoef', 'wbound'
                           'theta_z', 'theta_y', 'theta_x', 'tilt_y'])

    #d.update(ref_kw)
    try:
        ind = lines.index('*GRATINGNAME')
        grating_name = lines[ind + 1]
    except ValueError:
        grating_name = 'PRISM'

    #for line in lines:
    ind = lines.index('*THETAZ')
    d['theta_z'] = float(lines[ind + 1]) / 3600. # in degrees
    ind = lines.index('*THETAX')
    d['theta_x'] = float(lines[ind + 1]) / 3600. # in degrees
    ind = lines.index('*THETAY')
    d['theta_y'] = float(lines[ind + 1]) / 3600. # in degrees
    ind = lines.index('*TILTY')
    d['tilt_y'] = float(lines[ind + 1]) # in degrees
    try:
        ind = lines.index('*TILTX')
        d['tilt_x'] = float(lines[ind + 1]) # in degrees
    except ValueError:
        d['tilt_x'] = 0.0

    if disperser_type == 'gratingdata':
        ind = lines.index('*GROOVEDENSITY')
        d['groove_density'] = float(lines[ind + 1])
    elif disperser_type == 'prismdata':
        ind = lines.index('*ANGLE')
        d['angle'] = float(lines[ind + 1]) # in degrees

        ind = lines.index('*TREF')
        d['tref'] = float(lines[ind + 1]) # Temperature in K

        ind = lines.index('*PREF')
        d['pref'] = float(lines[ind + 1]) # Pressure in ATM
        for line in lines:
            if line.startswith('*COEFFORMULA'):
                ind = lines.index(line) + 1
                coefs = np.array(lines[ind : ind+int(line[-1])], dtype=np.float)
                kcoef = coefs[::2]
                lcoef = coefs[1::2]
                d['lcoef'] = list(lcoef)
                d['kcoef'] = list(kcoef)
            elif line.startswith('*THERMALCOEF'):
                # 6 coeffs - D0, D1, D2, E0, E1, lambdak
                ind = lines.index(line) + 1
                coefs = lines[ind : ind+int(line[-1])]
                d['tcoef'] = [float(c) for c in coefs]
            elif line.startswith('*WBOUND'):
                ind = lines.index(line) + 1
                coefs = [float(c) for c in (lines[ind]).split()]
                d['wbound'] = coefs

    assert grating_name in tiltyfile
    assert grating_name in tiltxfile
    tiltyd = disperser_tilt(tiltyfile)
    tiltxd = disperser_tilt(tiltxfile)

    d['gwa_tiltx'] = tiltyd
    d['gwa_tilty'] = tiltxd

    disperser_model = DisperserModel()
    if disperser_type == 'gratingdata':
        disperser_model.groovedensity = d['groove_density']
    else:
        disperser_model.angle = d['angle']
        disperser_model.kcoef = d['kcoef']
        disperser_model.lcoef = d['lcoef']
        disperser_model.tcoef = d['tcoef']
        disperser_model.pref = d['pref']
        disperser_model.tref = d['tref']

    disperser_model.gwa_tiltx = d['gwa_tiltx']
    disperser_model.gwa_tilty = d['gwa_tilty']
    disperser_model.theta_x = d['theta_x']
    disperser_model.theta_y = d['theta_y']
    disperser_model.theta_z = d['theta_z']
    disperser_model.tilt_x = d['tilt_x']
    disperser_model.tilt_y = d['tilt_y']
    disperser_model.wbound = d['wbound']
    disperser.meta.pgrating = "ANY|N/A|G140M|G140H|G235M|G235H|G395M|G395H|MIRROR|PRISM|"
    disperser_model.meta.exp_type = "N/A"



def disperser_tilt(tiltfile):
    with open(tiltfile) as f:
        s = f.read()
    tilt_d = {}
    vals = s.split('*')
    for line in vals[::-1]:
        if line.startswith("CoeffsTemperature00"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = {}
            idt_coef = [float(c) for c in l[1 : 1+n]]
            for i , c in enumerate(idt_coef[::-1]):
                coeffs['c' + str(i)] = c
            tilt_d['tilt_model'] = models.Polynomial1D(n-1, **coeffs)
        elif line.startswith("Temperatures"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = l[1:1+n]
            tilt_d['temperatures'] = [float(c) for c in coeffs]
        elif line.startswith("Zeroreadings"):
            l = line.split('\n')
            n = int(l[0].split()[1])
            coeffs = l[1:1+n]
            tilt_d['zeroreadings'] = [float(c) for c in coeffs]
        elif line.startswith("Unit"):
            tilt_d['unit'] = line.split('\n')[1]
    return tilt_d

if __name__ == '__main__':
    import argpars
    parser = argpars.ArgumentParser(description="Creates NIRSpec 'disperser' reference files in ASDF format.")
    parser.add_argument("disperser_dir", type=str, help="Directory with disperser files.")
    res = parser.parse_args()

    for i, grating in enumerate(["G140H", "G140M", "G235H", "G235M", "G395H",
                                 "G395M", "MIRROR", "PRISM"]):
        # disperser refers to the whole reference file, either on disk or asdf
        # grating to the parameter the reference file is selected on
        disperser_name =  "nirspec_disperser_{0}.asdf".format(grating)
        ref_kw = common_reference_file_keywords(reftype="disperser", title="NIRSPEC Disperser Description - CDP4",
                                                description="Grating as built description, tilt x/y/z fitted with FM2 CAL phase data.",
                                                exp_type="N/A", useafter=useafter, author=author, filename=disperser_name)
        ref_kw['grating'] = grating
        if grating == 'PRISM':
            dis_file = "disperser_" + grating + ".pri"
        else:
            dis_file = "disperser_" + grating + ".dis"
        gtpy_file = "disperser_" + grating + "_TiltY.gtp"
        gtpx_file = "disperser_" + grating + "_TiltX.gtp"
        disp_refname = os.path.join(ref_files, "Description", dis_file)
        tilty_refname = os.path.join(ref_files, "Description", gtpy_file)
        tiltx_refname = os.path.join(ref_files, "Description", gtpx_file)


        try:
            disperser2asdf(disp_refname, tilty_refname, tiltx_refname, disperser_name, ref_kw)
        except:
            raise Exception("Disperser file was not converted.")
