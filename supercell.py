from numpy import *
#from matplotlib.pyplot import *
from collections import *
from operator import itemgetter

vasp_data = []
vasp_file=raw_input("Write the name of input file(without .vasp):  ")
geo = raw_input("Is it free, surface, or periodic(f/s/p):  ")
if geo == 'f':
   geometry = 'free'
elif geo == 's':
   geometry = 'surface'
elif geo == 'p':
   geometry = 'periodic'
else:
   print "Type correct letter!"
   exit()

vasp_file_name = vasp_file +'.vasp'
with open(vasp_file_name, 'r') as position_data:
     for position_line in position_data:
          position_words = position_line.split()
          vasp_data.append(position_words)

ax, ay, az = vasp_data[2]
bx, by, bz = vasp_data[3]
cx, cy, cz = vasp_data[4]

lattice_a = sqrt(float(ax)**2 + float(ay)**2 + float(az)**2)
lattice_b = sqrt(float(bx)**2 + float(by)**2 + float(bz)**2)
lattice_c = sqrt(float(cx)**2 + float(cy)**2 + float(cz)**2)

types_of_atoms = vasp_data[5]
no_of_types_of_atoms = size(types_of_atoms)

no_of_atoms_in_each_type = vasp_data[6]
for i in range(0, no_of_types_of_atoms):
    no_of_atoms_in_each_type[i] = int(no_of_atoms_in_each_type[i])  # to convert int

xyz = []
xyz_data = []
x = []
y = []
z = []
no_of_atoms = sum(no_of_atoms_in_each_type)
for i in range(8,no_of_atoms+8):
    xyz = vasp_data[i]
    for j in range(0,3):
        xyz[j] = float(xyz[j])
    xyz_data.append(xyz)
xyz_data = array(xyz_data)
unit = ''.join(vasp_data[7])
x = xyz_data[:,0]
y = xyz_data[:,1]
z = xyz_data[:,2]

nx = input("Super cell value in x:",)
ny = input("Super cell value in y:",)
nz = input("Super cell value in z:",)

ax =float(ax) * nx
bx =float(bx) * ny
cx =float(cx) * nz
ay =float(ay) * nx
by =float(by) * ny
cy =float(cy) * nz
az =float(az) * nx
bz =float(bz) * ny
cz =float(cz) * nz

supercell_sy_xyz = []
no_of_atoms_in_supercell = nx*ny*nz*no_of_atoms
atomic_symbol = []
for i in range(0, no_of_types_of_atoms):
      for j in range(0,no_of_atoms_in_each_type[i]):
            atomic_symbol.append(types_of_atoms[i])
########################################################
             #Producing extra atoms for supercell
########################################################
for i in range(0, nx):
    for j in range(0, ny):
        for k in range(0, nz):
            for l in range(0, no_of_atoms):
                if unit in ('Direct', 'direct'):
                   supercell_sy_xyz.append([atomic_symbol[l], (x[l]+i)/nx, (y[l]+j)/ny, (z[l]+k)/nz])
                else:
                   supercell_sy_xyz.append([atomic_symbol[l], x[l]+lattice_a*i, y[l]+lattice_b*j, z[l]+lattice_c*k])
lattice_a = sqrt(float(ax)**2 + float(ay)**2 + float(az)**2)
lattice_b = sqrt(float(bx)**2 + float(by)**2 + float(bz)**2)
lattice_c = sqrt(float(cx)**2 + float(cy)**2 + float(cz)**2)
supercell_sy_xyz = sorted(supercell_sy_xyz, key=itemgetter(2))
atomic_symbol = []
x = []
y = []
z = []

for i in range (0, no_of_atoms_in_supercell):
   atomic_symbol.append(str(supercell_sy_xyz[i][0])) 
   x.append(float(supercell_sy_xyz[i][1]))
   y.append(float(supercell_sy_xyz[i][2]))
   z.append(float(supercell_sy_xyz[i][3]))

###############################################################
            # Printing output on terminal
###############################################################
types_of_atoms = list(set(atomic_symbol))
no_of_types_of_atoms = size(types_of_atoms)

no_of_atoms_in_each_type = []
for i in range (0, no_of_types_of_atoms):
    j = 0
    for k in range(0, no_of_atoms_in_supercell):
        if types_of_atoms[i] == atomic_symbol[k]:
           j += 1
    no_of_atoms_in_each_type.append(str(j))


print "Total number of atoms and unit :", no_of_atoms_in_supercell#, unit
print "Types of atoms                 :", ' '.join(types_of_atoms)
print "Number of atoms in each type   :", ' '.join(no_of_atoms_in_each_type)
print "Geometry & lattice             :", geometry, lattice_a, lattice_b, lattice_c

###############################################################
               # Printing vasp output
###############################################################
vacuum = raw_input("Do you want vacuum along z-axis? (Y/N):")
if vacuum == 'Y':
   if unit in ('Direct', 'direct'):
       print "Yes going"
       for i in range(0, no_of_atoms_in_supercell):
           z[i] = (z[i]*cz)/(cz+15+7.5)
   else:
       for i in range(0, no_of_atoms_in_supercell):
           z[i] = z[i] 
   cz = cz+15


outfile = raw_input("Write the name of output file that you want without .vasp:  ")
output_file2 = outfile+".vasp"
f1 = open(output_file2,"w")
print >>  f1, "supercell conversion by Selva"
print >> f1, "1"
print >> f1, "  ",(" %12.8f    %12.8f    %12.8f" %(ax, ay, az))
print >> f1, "  ",(" %12.8f    %12.8f    %12.8f" %(bx, by, bz))
print >> f1, "  ",(" %12.8f    %12.8f    %12.8f" %(cx, cy, cz))
print >> f1, ' '.join(types_of_atoms)
print >> f1, ' '.join(no_of_atoms_in_each_type)
print >> f1, unit
for i in range (0, no_of_types_of_atoms):
     no_of_atoms_in_each_type[i] = int(no_of_atoms_in_each_type[i])
for i in range(0, no_of_types_of_atoms):
        for k in range(0, no_of_atoms_in_supercell):
            if types_of_atoms[i] ==  atomic_symbol[k]:
                  print >> f1, (" %12.6f    %12.6f    %12.6f" %(x[k],y[k],z[k]))
###############################################################
               # Printing xyz output
###############################################################
#
#for i in range(0, no_of_atoms_in_supercell):
#    if atomic_symbol[i] == 'Cd':
#        atomic_symbol[i] = 'Cd_sc'
#
#for i in range(0, no_of_types_of_atoms):
#    if types_of_atoms[i] == 'Cd':
#        types_of_atoms[i] = 'Cd_sc'
#outfile1 = outfile+".xyz"
#f = open(outfile1,"w")
#print >> f, no_of_atoms_in_supercell, "angstroem"
#print >> f, geometry, lattice_a, lattice_b, lattice_c
#for k in range(0, no_of_atoms_in_supercell):
#     print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(atomic_symbol[k], x[k],  y[k],  z[k]))
###############################################################


#f.close()
f1.close()
exit()
