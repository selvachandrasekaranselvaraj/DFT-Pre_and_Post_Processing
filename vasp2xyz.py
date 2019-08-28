from numpy import *
from matplotlib.pyplot import *
vasp_data = []
vasp_file=raw_input("Write the name of input file(without .vasp):  ")
vasp_file_name = vasp_file +'.vasp'

geometry = raw_input("Is it free, surface, or periodic [defaut: periodic]:")
if geometry == '': geometry = 'periodic'

with open(vasp_file_name, 'r') as position_data:
     for position_line in position_data:
          position_words = position_line.split()
          vasp_data.append(position_words)

lattice_a, dummy, dummy = vasp_data[2]
dummy, lattice_b, dummy = vasp_data[3]
dummy, dummy, lattice_c = vasp_data[4]

lattice_a = float(lattice_a)
lattice_b = float(lattice_b)
lattice_c = float(lattice_c)

types_of_atoms = []
types_of_atoms = vasp_data[5]


no_of_atoms_in_each_type = []
no_of_atoms_in_each_type = vasp_data[6]
for i in range(0, size(vasp_data[6])):
    no_of_atoms_in_each_type[i] = int(no_of_atoms_in_each_type[i])
no_of_atoms = sum(no_of_atoms_in_each_type)

unit = str(vasp_data[7])
unit = unit[2:]
unit = unit[:-2]
xyz = []
xyz_data = []
x = []
y = []
z = []
for i in range(8,sum(no_of_atoms_in_each_type)+8):
    xyz = vasp_data[i]
    for j in range(0,3):
        xyz[j] = float(xyz[j])
    xyz_data.append(xyz)
xyz_data = array(xyz_data)
if unit in ('Direct', 'direct'):
   for i in range(no_of_atoms):
       if xyz_data[i,0] > 0.98:
          x.append(xyz_data[i,0]-1)
          y.append(xyz_data[i,1]*lattice_b)
          z.append(xyz_data[i,2]*lattice_c)
       elif xyz_data[i,1] > 0.98:
          x.append(xyz_data[i,0]*lattice_a)
          y.append(xyz_data[i,1]-1)
          z.append(xyz_data[i,2]*lattice_c)
       elif xyz_data[i,2] > 0.98:
          x.append(xyz_data[i,0]*lattice_a)
          y.append(xyz_data[i,1]*lattice_b)
          z.append(xyz_data[i,2]-1)
       else:
          x.append(xyz_data[i,0]*lattice_a)
          y.append(xyz_data[i,1]*lattice_b)
          z.append(xyz_data[i,2]*lattice_c)
else:
   for i in range(no_of_atoms):
       if xyz_data[i,0]/lattice_a > 0.98:
          x.append(xyz_data[i,0]-lattice_a)
          y.append(xyz_data[i,1])
          z.append(xyz_data[i,2])
       elif xyz_data[i,1]/lattice_b > 0.98:
          x.append(xyz_data[i,0])
          y.append(xyz_data[i,1]-lattice_b)
          z.append(xyz_data[i,2])
       elif xyz_data[i,2]/lattice_c > 0.98:
          x.append(xyz_data[i,0])
          y.append(xyz_data[i,1])
          z.append(xyz_data[i,2]-lattice_c)
       else:
          x.append(xyz_data[i,0])
          y.append(xyz_data[i,1])
          z.append(xyz_data[i,2])

#       x = xyz_data[:,0]
#       y = xyz_data[:,1]
#       z = xyz_data[:,2]

for i in range(0,sum(no_of_atoms_in_each_type)):
    if y[i] > lattice_b-1:
          y[i] = y[i] - lattice_b
    else:
          y[i] = y[i]

shift_y_z = "yes"
shift = raw_input("Do you want to shift y-z coordinates? (yes/no) [default: no]:")
if shift == '': shift = 'no'
outfile = vasp_file+'.xyz'
f = open(outfile,"w")
print >> f, sum(no_of_atoms_in_each_type), "angstroem"

atoms_to_include_semicore = ['Cd', 'Hg', 'Zn']
def semicore (atomic_symbol):
    if atomic_symbol in atoms_to_include_semicore:
       return atomic_symbol+'_sc'
    else:
       return atomic_symbol

no_of_types_of_atoms = size(types_of_atoms)
for i in range(0, no_of_types_of_atoms):
        types_of_atoms[i] = semicore(types_of_atoms[i])

k = 0
if shift_y_z == shift:
    print >> f, geometry, lattice_a, lattice_c, lattice_b
    for i in range(0, size(vasp_data[6])):
        for j in range(0,no_of_atoms_in_each_type[i]):
            print >> f, ("%-4s %11.6f    %11.6f    %11.6f" %(types_of_atoms[i], x[k],  z[k],  y[k]))
            k = k+1
    f.close()
else:
    print >> f, geometry, lattice_a, lattice_b, lattice_c
    for i in range(0, size(vasp_data[6])):
        for j in range(0,no_of_atoms_in_each_type[i]):
            print >> f, ("%-4s %11.6f    %11.6f    %11.6f" %(types_of_atoms[i], x[k],  y[k],  z[k]))
            k = k+1
    f.close()
