#!/usr/bin/env python2.7
"""
HullRad
Calculates hydrodynamic coefficients for a biomolecular structure.
AUTHOR: Patrick Fleming
DATE: 2017

If you publish work that uses this code please cite:
Fleming, PJ and Fleming, KG "HullRad: Fast Calculations of Folded and Disordered
Protein and Nucleic Acid Hydrogynamic Properties", Biophysical Journal, 2018
DOI: 10.1016/j.bpj.2018.01.002


The author would like to acknowledge the work of Jose Garcia de la Torre who led
the way in development of algorithms for calculating macromolecular
hydrodynamic properties.

--------------------------------------------------------------------------------
Requires python version2.7

Uses PDB file format as input.

Requires qconvex in path. This is a program in the qhull suite of programs.
Qhull may be downloaded from http://www.qhull.org/download.
	For UNIX download "Qhull_2015.2 for Unix"
	tar xzvf qhull-2015-src-7.2.0.tgz
	cd qhull-2015.2/
	make
	sudo make install
	and the path for qconvex will be /usr/local/bin/qconvex

OS X qhull binaries are available from MacPorts package manager. 
	type "sudo port install qhull".
	and the path for qconvex will be /opt/local/bin/qconvex

OS X binaries are also available from Fink package manager.

The default in this script is the UNIX path (/usr/local/bin/qconvex) but you may 
have to change it if you are using OS X.

After installing qhull go to: ### Edit the path to qconvex ### below (~line 476) and
change /usr/local/bin/qconvex to /opt/local/bin/qconvex if that is where qconvex is.

If you are using OS X but are familiar with UNIX and have gcc installed, the UNIX
instructions above work fine on a Mac.

FOR WINDOWS USERS:
There are precompiled versions of the qhull programs available for Windows 7 and up
at http://www.qhull.org/download
Unpack qhull in a directory of choice and ### Edit the path to qconvex ### below (~line 330)
to include the path where your executable was unpacked. 

Output looks like:
 
 pdb8rat.ent
 HYDROLASE (NUCLEIC ACID,RNA)            13-AUG-91   8RAT
  M               :        13692     g/mol
  v_bar           :        0.710     mL/g
  R(Anhydrous)    :        15.68     Angstroms
  Axial Ratio     :         1.31
  f/fo            :         1.21
  Dt              :       1.13e-06   cm^2/s
  R(Translation)  :        18.99     Angstroms
  s               :       1.84e-13   sec
  [eta]           :         3.16     cm^3/g
  Dr              :       1.91e+07   s^-1
  R(Rotation)     :        20.36     Angstroms

To write out the unified side chain psuedo-atom model uncomment the two lines
near the bottom labeled, # Write out model (~line 529)
"""

USAGE = """
python2.7 HullRad.py [PDBcode].pdb
    or 
./HullRad.py [PDBcode].pdb
"""

import os
import sys
import math
import string
import shlex, subprocess

# Amino Acid mass, volumes,partial specific volumes and unified side chain sphere radius
# Mass is minus water (from SEDNTERP database)
# Volumes from Cohn and Edsall (1943) as corrected by Perkins_EJB_1986 
# vbar is 0.60224(V/M)
#          Mass		Volume	   vbar		Radius of equiv. sphere of side chain volume
AA_data = {
 'ALA': (  71.0939,	87.2,      0.74,	1.852),
 'ARG': ( 156.203,	188.2,     0.73,	3.123),
 'ASN': ( 114.119,	120.1,     0.63,	2.422),
 'ASP': ( 115.104,	115.4,     0.60,	2.356),
 'CYS': ( 103.16,	106.7,     0.62,	2.224),
 'GLN': ( 128.16,	145.1,     0.68,	2.722),
 'GLU': ( 129.131,	140.9,     0.66,	2.676),
 'GLY': (  57.0669,	60.6,      0.64,	0.000),
 'HIS': ( 137.157,	152.4,     0.67,	2.798),
 'HSE': ( 137.157,	152.4,     0.67,	2.798),
 'HSP': ( 137.157,	152.4,     0.67,	2.798),
 'HSD': ( 137.157,	152.4,     0.67,	2.798),
 'ILE': ( 113.175,	168.9,     0.90,	2.957),
 'LEU': ( 113.175,	168.9,     0.90,	2.957),
 'LYS': ( 128.19,	174.3,     0.82,	3.005),
 'MET': ( 131.215,	163.1,     0.75,	2.903),
 'MSE': ( 178.125,	163.1,     0.75,	2.903),
 'PHE': ( 147.192,	187.9,     0.77,	3.121),
 'PRO': (  97.132,	122.4,     0.76,	2.453),
 'SER': (  87.093,	91.0,      0.63,	1.936),
 'THR': ( 101.12,	117.4,     0.70,	2.385),
 'TRP': ( 186.229,	228.5,     0.74,	3.422),
 'TYR': ( 163.192,	192.1,     0.71,	3.155),
 'VAL': (  99.148,	141.4,     0.86,	2.682)
}

# From voss_jmb_2005.pdf and nadassy_nar_2001.pdf
# Data for nucleoside, = nucleotide - 79.95 for PO3H
# Phosphate for na lacking terminal phosphates
# vbar is 0.60224(V/M)
#          Mass         Volume  vbar    Radius of equiv. sphere of base volume
NA_data = {
 '  A':  ( 267.25,       272.3,  0.614,   0.00),
 '1MA':  ( 281.25,       299.3,  0.641,   0.00),
 'MIA':  ( 383.45,       385.1,  0.605,   0.00),
 '  C':  ( 243.22,       248.1,  0.614,   0.00),
 '5MC':  ( 257.22,       275.1,  0.644,   0.00),
 'OMC':  ( 257.22,       275.1,  0.644,   0.00),
 '  G':  ( 283.24,       279.0,  0.593,   0.00),
 '2MG':  ( 297.27,       306.0,  0.620,   0.00),
 '7MG':  ( 299.27,       306.0,  0.620,   0.00),
 'M2G':  ( 311.27,       333.0,  0.644,   0.00),
 'OMG':  ( 297.27,       306.0,  0.620,   0.00),
 'YYG':  ( 508.54,       500.0,  0.593,   0.00),
 ' YG':  ( 428.55,       427.1,  0.600,   0.00),
 '  G':  ( 283.24,       279.0,  0.593,   0.00),
 '  U':  ( 244.20,       243.9,  0.602,   0.00),
 '5MU':  ( 258.20,       270.9,  0.632,   0.00),
 '4SU':  ( 260.35,       270.5,  0.626,   0.00),
 'H2U':  ( 246.20,       243.9,  0.602,   0.00),
 'PSU':  ( 244.20,       243.9,  0.602,   0.00),
 '  I':  ( 268.23,       273.0,  0.613,   0.00),
 ' DA':  ( 251.25,       271.0,  0.650,   0.00),
 ' DC':  ( 227.22,       246.8,  0.654,   0.00),
 ' DG':  ( 267.25,       277.7,  0.626,   0.00),
 ' DT':  ( 242.23,       264.4,  0.657,   0.00),
 ' DI':  ( 252.22,       271.7,  0.649,   0.00),
 'PO2':  (  62.97,	  52.4,  0.501,   0.00)
}

def distance(ax,ay,az, bx,by,bz):
    #Euclidean distance between two atoms
    return math.sqrt((ax - bx)**2.0 + (ay - by)**2.0 + (az - bz)**2.0)

def get_coords(model_array):
    # Get x,y,z coords of CA and CB atoms
    coords = []
    for row in range(len(model_array)):
        coords.append((float(model_array[row][4]), float(model_array[row][5]), float(model_array[row][6])))
    return coords

def model_from_pdb(file):
    # Makes a reduced atom  model of the pdb file
    # This is the model used to make the convex hull
    # For proteins the CB is displaced along the CA-CB vector a distance 
    #   equal to the radius of a sphere equal to the volume of the side chain
    # For nucleic acids only the oxygen and nitrogen atoms are used 
 
    # Read the PDB file
    data = open(file, 'rb').readlines()
 
    # Initialize
    # Initial protein atoms from PDB file
    prot_rec = []
    # Initial nucleic acid atoms from PDB file
    na_rec = []
    # Reduced list of atoms
    atom_rec = []
    struc_name = ' '
    num_MG = 0
    num_MN = 0

    # Get all relevant atoms even if in wrong order
    for line in data:
        # Get name of structure
        if (line[:6] == 'HEADER'):
            struc_name = line.strip()
            print struc_name[9:]

        elif (line[:4] == 'ATOM'):
            resname = line[17:20]
            if AA_data.has_key(resname) and (line[13:15] == 'N ' or \
              line[13:15] == 'CA' or \
              line[13:15] == 'C ' or \
              line[13:15] == 'O ' or \
              line[13:16] == 'OT1' or \
              line[13:15] == 'CB') and \
              (line[16] != 'B') and \
              (line[:3] != 'END'):
                prot_rec.append(line)
            elif NA_data.has_key(resname) and (line[13:14] == 'N' or \
              line[13:14] == 'O' or line[13:14] == 'P'):
                na_rec.append(line)
        elif (line[:6] == 'HETATM'):
            resname = line[17:20]
            if line[12:14] == 'MG':
                num_MG += 1
            elif line[12:14] == 'MN':
                num_MN += 1
            # Some nucleic acids are in HETATM (!)
            elif NA_data.has_key(resname) and (line[13:14] == 'N' or \
              line[13:14] == 'O' or line[13:14] == 'P'):
                na_rec.append(line)

    # Get only the CB atoms from the protein
    cb_rec = []
    for line in prot_rec:
        if (line[13:15] == 'CB') and \
           (line[0:3] != 'END'):
            cb_rec.append((line))

    # Count the CB atoms and add to atom list in correct order
    num_cb = len(cb_rec)
    line_counter = 0
    cb_counter = 0
    for line in prot_rec:
        if (line[13:15] == 'N ' or \
          line[13:15] == 'CA' or \
          line[13:15] == 'C ' or \
          line[13:15] == 'O ' or \
          line[13:16] == 'OT1'):
            atom_rec.append((line))
            line_counter += 1
            if (math.fmod(line_counter,4) == 0 and line[17:20] != 'GLY' \
                 and cb_counter < num_cb):
                atom_rec.append((cb_rec[cb_counter]))
                cb_counter +=1
            elif (math.fmod(line_counter,4) == 0 and line[17:20] != 'GLY' \
                 and cb_counter >= num_cb):
                print ' Found non-gly residue with no CB atom'
                sys.exit()

    # Put atom info into model array
    num_atoms =  len(atom_rec)
    model_array = [['X' for j in range(7)] for i in range(num_atoms)]
    for row in range(len(model_array)):
        model_array[row][0] = row
        model_array[row][1] = (atom_rec[row][11:16])
        model_array[row][2] = (atom_rec[row][17:20])
        model_array[row][3] = (atom_rec[row][22:26])
        model_array[row][4] = (atom_rec[row][30:38])
        model_array[row][5] = (atom_rec[row][38:46])
        model_array[row][6] = (atom_rec[row][46:54])

    # Make unified side chain as CB
    for row in range(len(model_array)):
        if (model_array[row][1] == '  CA ' and model_array[row][2] != 'GLY'):
            rnam = model_array[row][2]
            # CA coords
            ca_x = float(model_array[row][4])
            ca_y = float(model_array[row][5])
            ca_z = float(model_array[row][6])
            # CB coords (third next coord: C,O,CB)
            cb_x = float(model_array[row+3][4])
            cb_y = float(model_array[row+3][5])
            cb_z = float(model_array[row+3][6])
            # CA-CB distance
            dx = cb_x - ca_x
            dy = cb_y - ca_y
            dz = cb_z - ca_z
            ca_cb_dist = distance(ca_x,ca_y,ca_z,cb_x,cb_y,cb_z)
            # Extend CB outward along CA-CB vector
            aa_mass, vol_aa, vbar_aa, sc_rad_aa = AA_data[rnam]            
            new_cb_x = ca_x + (sc_rad_aa * dx/ca_cb_dist)
            new_cb_y = ca_y + (sc_rad_aa * dy/ca_cb_dist)
            new_cb_z = ca_z + (sc_rad_aa * dz/ca_cb_dist)
            model_array[row+3][4] = str(new_cb_x)
            model_array[row+3][5] = str(new_cb_y)
            model_array[row+3][6] = str(new_cb_z)

    # Put nucleic acid atoms into initial array
    num_atoms =  len(na_rec)
    na_array = [['X' for j in range(7)] for i in range(num_atoms)]
    for row in range(len(na_array)):
        na_array[row][0] = row
        na_array[row][1] = (na_rec[row][11:16])
        na_array[row][2] = (na_rec[row][17:20])
        na_array[row][3] = (na_rec[row][22:26])
        na_array[row][4] = (na_rec[row][30:38])
        na_array[row][5] = (na_rec[row][38:46])
        na_array[row][6] = (na_rec[row][46:54])

    # Add nucleic acid atoms to model array
    for row in range(len(na_array)):
        model_array.append(na_array[row])

    return num_MG,num_MN,model_array

def write_pdb(model_array, filename):
    # Write out reduced atom model in PDB format for display
    if type(filename) is type(''):
        ff = open(filename, 'wb')
    else:
        ff = filename

    write = ff.write

    pdbfmt = 'ATOM  %5d%3s %3s  %4d    %8.3f%8.3f%8.3f  0.00  0.00\n'
    for i in xrange(len(model_array)):
        a = model_array[i]
        write(pdbfmt % (int(a[0]),a[1],a[2],int(a[3]),float(a[4]),float(a[5]),float(a[6])))
    write('END\n')

    ff.close()

def Sved(num_MG,num_MN,model_array):
    #
    # Main function: Does most things and calls HullVolume
    #
    # Use Cohn-Edsall eq. to calc protein specific volume
    # vbar_prot = sum(ni*Mi*vbar_aa) / (ni*Mi) where n is number, M is mass
    #
    prot_mol_mass = 0.0
    numerator = 0.0
    NA = 0

    # Calc masses for protein and nucleic acid
    for row in range(len(model_array)):
        if (model_array[row][1] == '  CA '):
            rnam = model_array[row][2]
            aa_mass, vol_aa, vbar_aa, sc_rad_aa = AA_data[rnam]            
            numerator += (aa_mass * vbar_aa)
            prot_mol_mass += aa_mass
        elif ((model_array[row][1] == "  O4'") or (model_array[row][1] == "  O4*")):
            rnam = model_array[row][2]
            aa_mass, vol_aa, vbar_aa, sc_rad_aa = NA_data[rnam]
            numerator += (aa_mass * vbar_aa)
            prot_mol_mass += aa_mass
            NA = 1
        # Add PO2 because above were nucleosides
        elif (model_array[row][1] == "  P  "):
            numerator += (62.97 * 0.501)
            prot_mol_mass += 62.97

    # Add weight of water for N-term and C-term of protein
    prot_mol_mass += 18.0

    # If nucleic acid add weight of MG and MN
    ion_mass = (num_MG * 24.305) + (num_MN * 54.938)
    prot_mol_mass += ion_mass
     
    # Correct by -0.0025 because above volumes were measured at 25 deg C and want 20 deg C
    #  as per Svedberg, The Ultracentrifuge, 1940, Oxford Univ. Press
    vbar_prot = (numerator / prot_mol_mass) - 0.0025

    # Get convex hull area and volume 
    coords = get_coords(model_array)
    area_hull, vol_hull, Dmax = HullVolume(coords)

    # Calculate shell volume from variable of shell thickness
    # Hydration shell thickness is empirically determined to be optimal 
    # This expands each hull plane by hydration thickness
    vol_shell_wat = area_hull * 2.8

    # Calculate total volume of hydrated convex hull (including shell waters)
    vol_hyd_hull = vol_hull + vol_shell_wat

    # Estimate axial ratio 
    # Minus 3 Ang because DNA duplex a should be rod length, not diagonal
    #   of hull end vertices (Also necessary to make apoferritin axial ratio = 1)
    a = (Dmax/2.0) - 3.0
    b = math.sqrt((3.0 * vol_hull)/(4.0*math.pi*a))

    # But some very spherical structures may not have the diagonal a full 3 Ang longer.
    # Can't have a < b.
    if a > b: 
        # Translational Shape Factor
        numerator = math.sqrt(1.0 - (b/a)**2)
        denominator = ((b/a)**0.66666667) * math.log((1 + math.sqrt(1.0 - (b/a)**2))/(b/a))
        Ft = numerator/denominator
    else:
        a = b
        Ft = 1.0

    # Weight shape factor
    # Empirically found to work better with expanded volume
    #  (Many combinations tried)
    Ft = math.sqrt(Ft)

    # Axial ratio of prolate ellipsoid of same volume as convex hull
    a_b_ratio = a/b

    # Find radius of sphere of same volume as hydrated convex hull
    factor = 3.0/(4.0 * math.pi)
    Rht = (factor * vol_hyd_hull)**0.333333
    # Include Shape factor to give effective hydrodynamic translational radius
    Rht = Rht * Ft
    # Rht comes as Angstrom
    # need meters in equation below
    Rh_trans = Rht * 1e-8

    # Calculate Svedberg coeff.
    eta = 0.0100194 # poise
    rho = 0.998234 # g/ml water density
    fT = 6.0*math.pi*eta*Rh_trans
    s = prot_mol_mass * (1.0 - (vbar_prot * rho)) / (6.02214e23 * fT)

    # Calculate translational diffusion coeff.
    #R = 8.314e7
    kB = 1.381e-16
    T = 293.15
    Dt = kB*T / fT

    # Calculate fT/fo
    # Vbar is in ml/g, need A^3/Dalton
    vol_prot = prot_mol_mass*vbar_prot/0.60224
    Ro = ((3.0 * vol_prot) / (4.0 * math.pi))**0.333333
    ffo_hyd_P = Rht/Ro

    # Calculate rotational diffusion coeff.
    # Rotation more strongly affected by hydration, Halle & Davidovic, 2003
    # So increase hydration shell thickness
    # Hydration shell thickness empirically found to be optimal
    vol_shell_wat_rot = area_hull * 4.3
    vol_hyd_hull_rot = vol_hull + vol_shell_wat_rot
    Rhr = (factor * vol_hyd_hull_rot)**0.333333

    # Include empirical correction for shape factor to give effective 
    #  hydrodynamic rotational radius
    # This obtained by fitting to DNA duplexes
    Fr = Ft**4.0
    Rhr = Rhr * Fr
    # Rhr comes as Angstrom
    # need meters in equation below
    Rh_rot = Rhr * 1e-8
    fR = 8.0*math.pi*eta*(Rh_rot**3.0)
    Dr = kB*T / fR

    # Calculate intrinsic viscosity
    # From Einstein viscosity equation,
    # [n] = (2.5*6.02214e23*(4*pi*Rht^3)/3)/prot_mol_mass
    #
    # rearrange to
    # 
    int_vis = (10.0 * math.pi * 0.602214 * ((Rht)**3.0)) / (3.0 * prot_mol_mass)

    return s,Dt,Dr,vbar_prot,Rht,ffo_hyd_P,prot_mol_mass,Ro,Rhr,int_vis,a_b_ratio,Ft


def HullVolume(coords):
    #
    # Uses qconvex to calculate convex hull
    #
    #   coords are coords of reduced atom model
    #

    # Input for qconvex, etc
    dim = 3
    num_points = len(coords)
    instring = " %d\n " % (dim )
    instring += " %d\n " % (num_points )
    for i in range(len(coords)):
        instring += ' %f %f %f\n' % (coords[i][0],coords[i][1],
                                   coords[i][2])
    instring += '\n'

    # Call qconvex with argument, FA (compute total area and volume)
    ### Edit the path to qconvex ###
    p = subprocess.Popen(["/usr/local/bin/qconvex", "FA"],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    (input, output) = (p.stdin, p.stdout)
    input.write(instring)
    input.close()
    data = output.read()

    # Area of convex hull
    area_start = data.find("area:")
    area_end = data.find("Total volume:")
    if area_end == -1:
        area_end = data.find("Approximate volume:")
    area_hull = float(data[(area_start + 5) : area_end] )

    # Volume of convex hull
    vol_start = data.find("volume:")
    if vol_start == -1 :
        print data
        return 0
    vol_hull = float(data[vol_start + 7 :] )

    # Find max distance between vertices for Dmax calculation
    # Necessary for length of corresponding prolate ellipsoid of revolution
    # Run qconvex with different argument, p (returns vertex coordinates)
    p = subprocess.Popen(["/opt/local/bin/qconvex", "p"],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    (input, output) = (p.stdin, p.stdout)
    input.write(instring)
    input.close()
    data = map(string.split, output.readlines())
    Dmax = 0.0
    for i in range(2,len(data)):
        for j in range(3,len(data)):
            dist = math.sqrt((float(data[i][0])-float(data[j][0]))**2.0 + \
                             (float(data[i][1])-float(data[j][1]))**2.0 + \
                             (float(data[i][2])-float(data[j][2]))**2.0 )
            if dist > Dmax:
                Dmax = dist

    return area_hull, vol_hull, Dmax


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print USAGE
        sys.exit()

    file = sys.argv[1]

    # Print actual file name
    print(' %s' %  (file))

    # Make the reduce atom model
    num_MG,num_MN,model_array = model_from_pdb(file)

    # Write out model
    # Uncomment next two lines if you want the unified atom SC model written out
#   basename = os.path.splitext(file)[0]
#   write_pdb(model_array, basename + '_model.pdb')

    # Call the main function
    s,Dt,Dr,vbar_prot,Rht,ffo_hyd_P,M,Ro,Rhr,int_vis,a_b_ratio,Ft  = \
           Sved(num_MG,num_MN,model_array)

    # Print coeffiients to screen
    print('  M               :       %6.0f     g/mol' % (M))
    print('  v_bar           :       %6.3f     mL/g'  % (vbar_prot))
    print('  R(Anhydrous)    :       %6.2f     Angstroms' % (Ro))
    print('  Axial Ratio     :       %6.2f' % (a_b_ratio))
    print('  f/fo            :       %6.2f'  % (ffo_hyd_P))
    print('  Dt              :       %6.2e   cm^2/s' % (Dt))
    print('  R(Translation)  :       %6.2f     Angstroms' % (Rht))
    print('  s               :       %6.2e   sec' % (s))
    print('  [eta]           :       %6.2f     cm^3/g' % (int_vis))
    if a_b_ratio > 2.63:
        print(' ')
        print('  Caution. Axial ratio too large for accurate prediction')
        print('  of following rotational properties by HullRad.')

    print('  Dr              :       %6.2e   s^-1' % (Dr))
    print('  R(Rotation)     :       %6.2f     Angstroms' % (Rhr))
