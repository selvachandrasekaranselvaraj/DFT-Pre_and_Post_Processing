from numpy import *
from matplotlib.pyplot import *
from collections import *
from operator import itemgetter
import os


def change_coordinates(lat_a, lat_b, lat_c, x_f, y_f, z_f, lat_a_new):
  x_new=(x_f/lat_a*lat_a_new)
  y_new=(y_f/lat_b*lat_a_new)
  z_new=(z_f/lat_c*lat_a_new)
  return x_new, y_new, z_new
def change_lattice(lat_a, lat_b, lat_c, lat_a_new):
    lat_a = lat_a_new
    lat_b = lat_a_new
    lat_c = lat_a_new
    return lat_a, lat_b, lat_c

def shiftx(a):
    for i in range(no_of_atoms):
        x[i] += a
    return x


xyz_data = []
xyz_file=raw_input("Write the name of input file without .xyz  : ")
xyz_file_name=xyz_file+".xyz"
with open(xyz_file_name, 'r') as position_data:
     for position_line in position_data:
          position_words = position_line.split()
          xyz_data.append(position_words)

line_1,line_2 = xyz_data[0:2]
no_of_atoms = int(line_1[0])
unit        = str(line_1[1])
geometry   = str(line_2[0])
A = "atomic"
B = "atomicd0"
C = "Bohr"
if unit in ("atomic", "atomicd0", "Bohr"):
    lattice_a   = 0.529177*float(line_2[1])
    lattice_b   = 0.529177*float(line_2[2])
    lattice_c   = 0.529177*float(line_2[3])
else:
    lattice_a   = float(line_2[1])
    lattice_b   = float(line_2[2])
    lattice_c   = float(line_2[3])

atomic_symbol = []
x = []
y = []
z = []
for i in range(0,no_of_atoms):
    dummy = xyz_data[i+2]
    atomic_symbol.append(str(dummy[0]))
    if unit in ("atomic", "atomicd0", "Bohr"):
        x.append(0.529177*float(dummy[1]))
        y.append(0.529177*float(dummy[2]))
        z.append(0.529177*float(dummy[3]))
    else:
        x.append(float(dummy[1]))
        y.append(float(dummy[2]))
        z.append(float(dummy[3]))

types_of_atoms = list(set(atomic_symbol))
no_of_types_of_atoms = size(types_of_atoms)

no_of_atoms_in_each_type = []
for i in range (0, no_of_types_of_atoms):
    j = 0
    for k in range(0, no_of_atoms):
        if types_of_atoms[i] == atomic_symbol[k]:
           j += 1
    no_of_atoms_in_each_type.append(str(j))


atoms_to_include_semicore = ['Cd', 'Hg', 'Zn']
def semicore (atomic_symbol):
    if atomic_symbol in atoms_to_include_semicore:
       return atomic_symbol+'_sc'
    else:
       return atomic_symbol

for i in range(0, no_of_atoms):
        atomic_symbol[i] = semicore(atomic_symbol[i])

for i in range(0, no_of_types_of_atoms):
        types_of_atoms[i] = semicore(types_of_atoms[i])

###############################################################
            # Printing output on terminal
###############################################################
print
print "Total number of atoms and unit :", no_of_atoms, unit
print "Types of atoms                 :", ' '.join(types_of_atoms)
print "Number of atoms in each type   :", ' '.join(no_of_atoms_in_each_type)
print "Geometry & lattice             :",geometry, lattice_a, lattice_b, lattice_c
print
changes = input("Do you want to change LATTICE CONSTANT(1), ATOMIC TYPE(2), BOTH(3), or shift_coordinates(4) (1/2/3/4): ",)
if changes == 1:
   print
   lattice_a_new = input("Give your new lattice value: ",)
   for i in range (0, no_of_atoms):
       x[i], y[i], z[i]  = change_coordinates(lattice_a, lattice_b, lattice_c, x[i], y[i], z[i], lattice_a_new)

   lattice_a, lattice_b, lattice_c = change_lattice(lattice_a, lattice_b, lattice_c, lattice_a_new)
   print "Geometry & lattice             :",geometry, lattice_a, lattice_b, lattice_c

elif changes == 2:
   for i in range(0, no_of_types_of_atoms):
       print "Do you want to change atom", "\""+types_of_atoms[i]+"\"", "(y/n): "
       atom_change = raw_input()
       if atom_change == 'y':
          types_of_atoms_new = raw_input("Give new atomic symbol: ",)
          for j in range(0, no_of_atoms):
             if atomic_symbol[j] == types_of_atoms[i]:
                atomic_symbol[j] = types_of_atoms_new
       print

elif changes == 3:
   print
   lattice_a_new = input("Give your new lattice value: ",)
   for i in range (0, no_of_atoms):
       x[i], y[i], z[i]  = change_coordinates(lattice_a, lattice_b, lattice_c, x[i], y[i], z[i], lattice_a_new)

   lattice_a, lattice_b, lattice_c = change_lattice(lattice_a, lattice_b, lattice_c, lattice_a_new)
   print "Geometry & lattice             :",geometry, lattice_a, lattice_b, lattice_c
   for i in range(0, no_of_types_of_atoms):
       print "Do you want to change atom", "\""+types_of_atoms[i]+"\"", "(y/n): "
       atom_change = raw_input()
       if atom_change == 'y':
          types_of_atoms_new = raw_input("Give new atomic symbol: ",)
          for j in range(0, no_of_atoms):
             if atomic_symbol[j] == types_of_atoms[i]:
                atomic_symbol[j] = types_of_atoms_new
       print
elif changes == 4:
    print
    make_atom_center = input("Do you want to shift an atom to center(1), origin(2), or else(3) : (1/2/3)", )
    if make_atom_center == 1:
        atom_no = input("Atomic number:", )
        x_shift, y_shift, z_shift = x[atom_no-1], y[atom_no-1], z[atom_no-1]
        print "shifting atom: ", x_shift, y_shift, z_shift
        x_shift, y_shift, z_shift = (lattice_a/2.)-x[atom_no-1], (lattice_a/2.)-y[atom_no-1], (lattice_a/2.)-z[atom_no-1]
        print "shifting atom: ", x_shift, y_shift, z_shift
        for i in range(no_of_atoms):
            x[i] += x_shift
            y[i] += y_shift
            z[i] += z_shift
    elif make_atom_center == 2:
        atom_no = input("Atomic number:", )
        x_shift, y_shift, z_shift = x[atom_no-1], y[atom_no-1], z[atom_no-1]
        print "shift atom at: ", x_shift, y_shift, z_shift
        for i in range(no_of_atoms):
            x[i] += x_shift
            y[i] += y_shift
            z[i] += z_shift
    elif make_atom_center == 3:
        shift = raw_input("Do you want to shift along x, y, z, xy, xz, zx, xyz: ", )
        if shift == 'x':
            x_shift = input("Give the shift value along x: ", )
            x = shiftx(x_shift)
        elif shift == 'y':
            x_shift = input("Give the shift value along y: ", )
            y = shiftx(x_shift)
        elif shift == 'z':
            x_shift = input("Give the shift value along z: ", )
            z = shiftx(x_shift)
###############################################################
               # Printing output
###############################################################
cmd = "cp "+str(xyz_file)+".xyz "+str(xyz_file)+"_old.xyz"
os.system(cmd)
cmd = "Write the name of output file that you want without .xyz[default: "+str(xyz_file)+"]:  "
outfile = raw_input(cmd)
if outfile == '': outfile = str(xyz_file)+".xyz"
else: outfile = outfile+".xyz"

def align_coordinate(coordinate, lattice):
    if coordinate > lattice - 0.3:
       coordinate -= lattice
       return coordinate
    else:
       return coordinate

f = open(outfile,"w")
print >> f, no_of_atoms, "angstroem" #, energy, energy_unit
print >> f, geometry, lattice_a, lattice_b, lattice_c
if geometry == 'surface' or geometry == 'free':
  for k in range(0, no_of_atoms):
     x[k] = align_coordinate(x[k], lattice_a)
     z[k] = align_coordinate(z[k], lattice_c)
     print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(atomic_symbol[k], x[k],  y[k],  z[k]))
else:
  for k in range(0, no_of_atoms):
     x[k] = align_coordinate(x[k], lattice_a)
     y[k] = align_coordinate(y[k], lattice_b)
     z[k] = align_coordinate(z[k], lattice_c)
     print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(atomic_symbol[k], x[k],  y[k],  z[k]))
f.close()
