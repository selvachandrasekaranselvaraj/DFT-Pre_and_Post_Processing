from numpy import *
from matplotlib.pyplot import *
from collections import *
from operator import itemgetter, attrgetter

xyz_data = []
file_name=raw_input("Write the name of input file without .xyz  : ")
xyz_file_name = file_name+".xyz"
with open(xyz_file_name, 'r') as position_data:
     for position_line in position_data:
          position_words = position_line.split()
          xyz_data.append(position_words)

line_1,line_2 = xyz_data[0:2]
no_of_atoms = int(line_1[0])
unit        = str(line_1[1])
geometry   = str(line_2[0])

if unit in ["atomic", "atomicd0", "Bohr"]:
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
    if unit in ["atomic", "atomicd0", "Bohr"]:
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


print "Total number of atoms and unit :", no_of_atoms, unit
print "Types of atoms                 :", ' '.join(types_of_atoms)
print "Number of atoms in each type   :", ' '.join(no_of_atoms_in_each_type)
print "Geometry & lattice             :",geometry, lattice_a, lattice_b, lattice_c


#output_file = raw_input("Give name of of output file without .vasp   : ")
output_file = file_name+ ".vasp"
f = open(output_file,"w")
print >>  f, "xyz to vasp conversion by Selva"
print >> f, "1"
dummy = 0.0000000
print >> f, "  ",(" %12.8f    %12.8f    %12.8f" %(lattice_a, dummy, dummy))
print >> f, "  ",(" %12.8f    %12.8f    %12.8f" %(dummy, lattice_b, dummy))
print >> f, "  ",(" %12.8f    %12.8f    %12.8f" %(dummy, dummy, lattice_c))
print >> f, ' '.join(types_of_atoms)
print >> f, ' '.join(no_of_atoms_in_each_type)
print >> f, "Cartesian"


for i in range(0, no_of_types_of_atoms):
        for k in range(0, no_of_atoms):
            if types_of_atoms[i] ==  atomic_symbol[k]:
                  print >> f, (" %12.6f    %12.6f    %12.6f" %(x[k],y[k],z[k]))

f.close()
