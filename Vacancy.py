from numpy import *
#from matplotlib.pyplot import *
from collections import *
from operator import itemgetter

def pristine_bond_length(no_of_types_of_atoms, types_of_atoms):
    atom_pairs = ["C-C", "Si-C", "Si-Si", "Si-Ge", "Ge-Ge", "Cd-S", "Cd-Te", "Ga-As", "Ge-C", "Se-Te", "Se-Cd"]
    bulk_bond_lengths = [1.529, 1.868, 2.33117, 2.373, 2.416, 2.497, 2.778, 2.396, 1.948, 2.8, 2.8]
    if no_of_types_of_atoms >1:
       atom_pair = str(types_of_atoms[0]+'-'+types_of_atoms[1])
       if atom_pair in atom_pairs:
          return bulk_bond_lengths[atom_pairs.index(atom_pair)]
       else:
          atom_pair = str(types_of_atoms[1]+'-'+types_of_atoms[0])
          if atom_pair in atom_pairs:
             return bulk_bond_lengths[atom_pairs.index(atom_pair)]
          else:
             print "This atomic_pair is not available"
    else:
       try:
          atom_pair = str(types_of_atoms[0]+'-'+types_of_atoms[0])
          if atom_pair in atom_pairs:
             return bulk_bond_lengths[atom_pairs.index(atom_pair)]
       except:
         print "This atomic_pair is not available"


atoms_to_include_semicore = ['Cd', 'Hg', 'Zn']
def semicore (atomic_symbol):
    if atomic_symbol in atoms_to_include_semicore:
       return atomic_symbol+'_sc'
    else:
       return atomic_symbol

def distance_btw_atoms(A,B):
    x1, y1, z1 = A[1], A[2], A[3]
    x2, y2, z2 = B[1], B[2], B[3]
    distance =  sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
    return round(distance, 3)




xyz_data = []
xyz_file= raw_input("Write the name of input file without .xyz  : ")
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
    symbol = (dummy[0])
    if len(symbol) > 2:
        symbol = symbol[0:2]
    atomic_symbol.append(symbol)
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
bulk_bond_length = pristine_bond_length(no_of_types_of_atoms, types_of_atoms)

no_of_atoms_in_each_type = []
for i in range (0, no_of_types_of_atoms):
    j = 0
    for k in range(0, no_of_atoms):
        if types_of_atoms[i] == atomic_symbol[k]:
           j += 1
    no_of_atoms_in_each_type.append(int(j))
for i in range(0, no_of_atoms):
        atomic_symbol[i] = semicore(atomic_symbol[i])

print
work = raw_input("Do you want to ANALYSIS or CREATE vacancy (A/C)[default: A]:" )
if work == '': work = 'A'

###############################################################
#Sorting y-coordinates
sy_xyz = []
for i in range (0,no_of_atoms):
    dummy = [atomic_symbol[i], x[i], y[i], z[i]]
    sy_xyz.append(dummy)
sy_xyz = sorted(sy_xyz, key=itemgetter(2))

vacancy_atom = ['Xx', lattice_a/2, lattice_a/2, lattice_a/2]
tetrahedron_edge = lattice_a*sqrt(2)/6
###############################################################
print bulk_bond_length
fnn = []  # First nearest neighbor from vacancy
fnn_index = []
fnn_dis = []
snn = []
snn_index = []
snn_dis = []
nn_index = []
for j in range(0,no_of_atoms):
    dis = distance_btw_atoms(sy_xyz[j], vacancy_atom) 
    dis = (round(dis,5))
    if dis < bulk_bond_length*1.1:  
       dummy = sy_xyz[j]
       fnn.append(dummy)
       fnn_dis.append(dis)
       fnn_index.append(j)
       nn_index.append(j)
    elif dis > bulk_bond_length + 0.2 and dis < tetrahedron_edge + 0.2:
       dummy = sy_xyz[j]
       snn.append(dummy)
       snn_index.append(j)
       snn_dis.append(dis)
       nn_index.append(j)
print len(fnn), len(snn)
print "No. of FNN & SNN: ", len(fnn), len(snn)

fnn_index.sort(reverse=True)
snn_index.sort(reverse=True)
nn_index.sort(reverse=True)
for i in range(len(nn_index)):
    del(sy_xyz[nn_index[i]]) 
for i in range(len(fnn_index)):
    if fnn_dis[i] < 0.5:
       del(fnn[i])
       del(fnn_index[i])

###############################################################
dx = []  # FNN distances
dx.append(distance_btw_atoms(fnn[0], fnn[1]))
dx.append(distance_btw_atoms(fnn[0], fnn[2]))
dx.append(distance_btw_atoms(fnn[0], fnn[3]))
dx.append(distance_btw_atoms(fnn[1], fnn[2]))
dx.append(distance_btw_atoms(fnn[1], fnn[3]))
dx.append(distance_btw_atoms(fnn[2], fnn[3]))
dx = sorted(dx)
te = tetrahedron_edge
print "Pristine Edge: ", ("%2.3f" %(te))
print "Relaxed Edges: ", ("%2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f" %(dx[0], dx[1], dx[2], dx[3], dx[4], dx[5]))
print "Ratio        : ", ("%2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f" %( dx[0]/te, dx[1]/te, dx[2]/te, dx[3]/te, dx[4]/te, dx[5]/te))
#print ("%2.3f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f" %(te, dx[0]/te, dx[1]/te, dx[2]/te, dx[3]/te, dx[4]/te, dx[5]/te))
print
###############################################################


if work == 'A':
   exit()
print
print "Possible vaancy symmetries are Td_up(1), Td_down(2), D2d(3), C2v(4), C3v(5)"
symmetry = input("Type vacancy symmetry(1/2/3/4/5): ", )

def xyz_outward_move(A, B, r):  # Move coordinates toward vacancy
   C = []
   C.append(B[0])
   C.append(B[1] - (A[1]-B[1])/r)
   C.append(B[2] - (A[2]-B[2])/r)
   C.append(B[3] - (A[3]-B[3])/r)
   return C


def xyz_inward_move(A, B, r):  # Move coordinates toward vacancy
   C = []
   C.append(B[0])
   C.append(B[1] + (A[1]-B[1])/r)
   C.append(B[2] + (A[2]-B[2])/r)
   C.append(B[3] + (A[3]-B[3])/r)
   return C

def xy_inward_move(A, B, r):  # Move coordinates toward vacancy
   C = []
   C.append(B[0])
   C.append(B[1] + (A[1]-B[1])/r)
   C.append(B[2] + (A[2]-B[2])/r)
   C.append(B[3])
   return C

def Td_down():
    ratio = 100
    d_fnn = te #
    while d_fnn > bulk_bond_length + (bulk_bond_length*5/100):
          fnn = []
          for i in range(len(fnn_original)):   
              fnn.append(xyz_inward_move(vacancy_atom, fnn_original[i], ratio))
          d_fnn = distance_btw_atoms(fnn[0], fnn[1])
          ratio -= 0.01
    snn = []
    for i in range(len(snn_original)):
        snn.append(xyz_inward_move(vacancy_atom, snn_original[i], ratio*4))
    return fnn + snn


def Td_up():
    fnn = []
    for i in range(len(fnn_original)):
        fnn.append(xyz_outward_move(vacancy_atom, fnn_original[i], 5))
    snn = []
    for i in range(len(snn_original)):
        snn.append(xyz_outward_move(vacancy_atom, snn_original[i], 20))
    return fnn + snn


def D2d():
    ratio = 100
    d_fnn = te #
    while d_fnn > bulk_bond_length:
          fnn = []
          for i in range(len(fnn_original)):   
              fnn.append(xy_inward_move(vacancy_atom, fnn_original[i], ratio))
          d_fnn_t = []    
          for i in range(4):
              for j in range(4):
                  if i != j and i < j:
                     d_fnn_t.append(distance_btw_atoms(fnn[i], fnn[j]))
          d_fnn = min(d_fnn_t)
          ratio -= 0.01
    snn = []
    for i in range(len(snn_original)):
        snn.append(xy_inward_move(vacancy_atom, snn_original[i], ratio*4))
    return fnn + snn


def C2v():
    ratio = 100
    d_fnn = te #
    while d_fnn > bulk_bond_length:
          fnn = []
          for i in range(4):
              if vacancy_atom[3] > fnn_original[i][3]:
                  fnn.append(xy_inward_move(vacancy_atom, fnn_original[i], ratio))
              else:
                  fnn.append(xyz_outward_move(vacancy_atom, fnn_original[i], 5))
          d_fnn_t = []
          for i in range(4):
              for j in range(4):
                  if i != j and i < j:
                     d_fnn_t.append(distance_btw_atoms(fnn[i], fnn[j]))
          d_fnn = min(d_fnn_t)
          ratio -= 0.01
    snn = []
    for i in range(len(snn_original)):
        if vacancy_atom[3] > snn_original[i][3] and snn_original[i][2] == snn_original[i][1]:
            snn.append(xy_inward_move(vacancy_atom, snn_original[i], ratio*4))
        else:
            snn.append(xyz_outward_move(vacancy_atom, snn_original[i], 20))
    return fnn + snn


def C3v():
    d_fnn = te #
    fnn = []
    fnn.append(xyz_outward_move(vacancy_atom, fnn_original[0], 6))
    for i in range(1, 4):
        fnn.append(xyz_inward_move(vacancy_atom, fnn_original[i], 30))
    snn = []
    for i in range(len(snn_original)):
        dis = distance_btw_atoms(fnn_original[0], snn_original[i])
        if dis < bulk_bond_length:
            snn.append(xyz_outward_move(vacancy_atom, snn_original[i], 7))
        else:
            snn.append(xyz_outward_move(vacancy_atom, snn_original[i], 40))
    return fnn + snn

fnn_original = fnn
snn_original = snn
if symmetry == 1:
    nn = Td_up()
    symmetry = 'Td_up'
elif symmetry == 2:
    nn = Td_down()
    symmetry = 'Td_down'
elif symmetry == 3:
    nn = D2d()
    symmetry = 'D2d'
elif symmetry == 4:
    nn = C2v()
    symmetry = 'C2v'
elif symmetry == 5:
    nn = C3v()
    symmetry = 'C3v'
else:
    print "type correct symmetry"
    exit()


outfile = str(symmetry)+"_"+xyz_file+".xyz"
f = open(outfile,"w")
no_of_atoms = len(sy_xyz) + len(nn)
print >> f, no_of_atoms, "angstroem"
print >> f, geometry, lattice_a, lattice_b, lattice_c
for i in range(len(nn)):
    print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(nn[i][0], nn[i][1], nn[i][2], nn[i][3])) 


for k in range(len(sy_xyz)):
    print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(sy_xyz[k][0], sy_xyz[k][1], sy_xyz[k][2], sy_xyz[k][3]))
print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %('Se', vacancy_atom[1], vacancy_atom[2], vacancy_atom[3])), "# vacancy atom"
f.close()
