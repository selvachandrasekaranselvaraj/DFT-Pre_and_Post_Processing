from operator import itemgetter, attrgetter
from numpy import *
from matplotlib.pyplot import *
from collections import *
import re

import os
import commands
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import math

#print "Check the no. of node you want in pbs file!"
no_of_nodes = 2 #input("How many No_of_nodes are needed in pbs file?: ")
#Q1 = raw_input("Do you want to check now? (Y/N): ")
#if Q1 == 'Y':
#   exit()
print "******************************************************* \n"

################ Reading yaml file #############################
input_file = 'default.yaml'
try:
   inp_data = open(input_file, 'r')
except:
   print "There is no default.yaml file"
   exit()
for fline in inp_data:
    if 'nspin' in fline.split():
       dummy = fline.split()
       nspin = int(dummy[2])
    if 'mpol' in fline.split():
       dummy = fline.split()
       mpol = dummy[2]
    if 'ixc' in fline.split():
       dummy = fline.split()
       ixc  = dummy[2] 
       if ixc in ['1', 'LDA']:
          print "You are using LDA FUNCTIONAL"
       elif ixc in ['11', 'PBE']:
          print "You are using PBE FUNCTIONAL"

print "**********************************************************\n"
################ Reading xyz file #############################
#print "Searching xyz files..."
#xyz_file_names = []
#cmd1 = "ls *xyz"
#status, output = commands.getstatusoutput(cmd1)
#xyz_file_name = output.split()
#dummy = xyz_file_name[0] 
#dummy = dummy[-3:]
#if dummy == 'xyz':
#  for i in range(0, size(xyz_file_name)):
#     dummy1 = xyz_file_name[i]
#     xyz_file_names.append((dummy1))
#  xyz_file_names = xyz_file_names
#else:
#  print "There is no .xyz file in the directory"
#  exit()
#print "No. of xyz files found are : ", size(xyz_file_names)
#print "The names of xyz files found are : ", ', '.join(xyz_file_names)
#print "**********************************************************\n"
#################################################################
          #Reading xyz_data 
#################################################################
xyz_data = []
#xyz_file_name = raw_input("Write the name of input xyz file that you wnat to use (without .xyz): ")
xyz_file_name = 'opt'
xyz_file_name  = xyz_file_name+'.xyz'
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
def semicore (atomic_symbol, ixc):
    if ixc == 'LDA' or ixc == 1:
       if atomic_symbol in atoms_to_include_semicore:
          return atomic_symbol+'_sc'
       else:
          return atomic_symbol

    else:
       return atomic_symbol[0:2]

for i in range(0, no_of_atoms):
        atomic_symbol[i] = semicore(atomic_symbol[i], ixc)

for i in range(0, no_of_types_of_atoms):
        types_of_atoms[i] = semicore(types_of_atoms[i], ixc)


#print "******************************************************* \n"
#print " Name of types of atoms        : ", ', '.join(set(atomic_symbol)) 
#print " No. of atoms in each type     : ", ', '.join(no_of_atoms_in_each_type)
#print " Total No. of atoms            : ", no_of_atoms
lattice = [str(lattice_a), str(lattice_b), str(lattice_c)]
#print " Lttice constants a, b, c are  : ", ', '.join(lattice)




#################################################################
            # Writing Function for printing xyz output
#################################################################
if nspin == 2:
   print "*******************************************************\n"
   print "You are calculating spin polarized calculation with mpol =", mpol
   fsymbol1 = raw_input("Therefore, Write atomic symbols that you want to apply spin: ")
   print "******************************************************* \n"
   fsymbol2 = fsymbol1.split()
else:
   fsymbol2 = 'Xx'
   
def xyz_printing(fsymbol, fx, fy, fz, fmm):
   if fsymbol in fsymbol2:
       print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(fsymbol, fx, fy, fz)), fmm
   else:
       print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(fsymbol, fx, fy, fz)) 
#################################################################
            # Writing Function for printing pbs files
#################################################################
job_name = raw_input("Give job name to inclulde in pbs file: ")
def pbs_printing(f_xyz_file_name, f_pbs_file_names):
       f_xyz_file_name = outfile[:-4]
       f_pbs_file = f_xyz_file_name+".pbs"
       f_pbs_file_names.append(f_pbs_file)
       f2 = open(f_pbs_file, "w")
       print >> f2, "#!/bin/bash -l"
       print >> f2, "#PBS -q csp"
       print >> f2, "#PBS -N", job_name+'-'+f_xyz_file_name
       print >> f2, "#PBS -j oe"
       #nodes = floor(no_of_processors/32)
       #nodes = int(nodes)+1
       nodes = str(no_of_nodes)
       cmd3 = '#PBS -l nodes='+nodes+':ppn=32'
       print >> f2, cmd3
       #print >> f, "#PBS -l nodes=",nodes,":ppn=32"
       print >> f2, "#PBS -l walltime=24:00:00"
       print >> f2, "##PBS -l pmem=2000mb"
       print >> f2, "module load icc ompi mkl"
       print >> f2, "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/W/dc200827/bigdft-suite/tmp/install/lib"
       print >> f2, "BIGDFT=/W/dc200827/bigdft-suite/tmp/install/bin/bigdft"
       print >> f2, "export OMP_NUM_THREADS=1"
       print >> f2, "cd $PBS_O_WORKDIR"
       print >> f2, "mpirun $BIGDFT -l no -n",f_xyz_file_name,">",f_xyz_file_name+".log"
       f2.close()

###################################################################
apply_strain = 'Y' #raw_input("Do you want apply strain? (Y/N): ")
print "******************************************************* \n"
pbs_file_names = []
if apply_strain == 'N':
###################################################################
               # Output xyz file for without strain
###################################################################
    outfile = "".join(set(atomic_symbol))+'.xyz'
    f = open(outfile,"w")
    print >> f, no_of_atoms, "angstroem"
    print >> f, geometry, lattice_a, lattice_b, lattice_c
    for j in range(0,no_of_atoms):
        xyz_printing(atomic_symbol[j], x[j],  y[j],  z[j], mpol)
    pbs_printing(xyz_file_name, pbs_file_names)
    f.close( )

#####################################################################
elif apply_strain == 'Y':
#####################################################################
           # Getting informations  for applying strain  
#####################################################################
    print "E : Appling strain equally to axes"
    print "P : Appling strian like permutation"
    strain_axis = raw_input("Do you want apply strain along 'x' or 'y' or 'z' or 'xy-E' or 'xy-P' or 'yz' or 'zx' or 'xyz-E' or 'xyz-P''? : ")
    print "******************************************************* \n"
    def strained_lattice(lattice_I, strain_i):
        return lattice_I * (1.0 + float(strain_i) / 100.0)

    def strained_coordinate(coordinate_I, strain_i):
        return coordinate_I * (1.0 + (float(strain_i) / 100.0))

##########################  Apply strain  #########################
    if strain_axis == 'x':
       for i in range(-4, 5):
           i = round(float(i)/1.0, 1)
           outfile = str(i)+'.xyz'
           f = open(outfile,"w")
           print >> f, no_of_atoms, "angstroem"
           print >> f, geometry, strained_lattice(lattice_a, i), lattice_b, lattice_c
           for j in range(0,no_of_atoms):
               xyz_printing(atomic_symbol[j], strained_coordinate(x[j], i),  y[j],  z[j], mpol)
           pbs_printing(xyz_file_name, pbs_file_names)
       f.close( )

     
    elif strain_axis == 'y':
       for i in range(-4, 5):
           i = round(float(i)/1.0, 1)
           outfile = str(i)+'.xyz'
           f = open(outfile,"w")
           print >> f, no_of_atoms, "angstroem"
           print >> f, geometry, lattice_a, strained_lattice(lattice_b, i), lattice_c
           for j in range(0,no_of_atoms):
               xyz_printing(atomic_symbol[j], x[j], strained_coordinate(y[j], i),  z[j], mpol)
           pbs_printing(xyz_file_name, pbs_file_names)
       f.close( )


    elif strain_axis == 'z':
       for i in range(-4, 5):
           i = round(float(i)/1.0, 1)
           outfile = str(i)+'.xyz'
           f = open(outfile,"w")
           print >> f, no_of_atoms, "angstroem"
           print >> f, geometry, lattice_a, lattice_b, strained_lattice(lattice_c, i)
           for j in range(0,no_of_atoms):
               xyz_printing(atomic_symbol[j], x[j], y[j], strained_coordinate(z[j], i), mpol)
           pbs_printing(xyz_file_name, pbs_file_names)
       f.close( )
     
     
    elif strain_axis == 'xy-E':
       for i in range(-4, 5):
           i = round(float(i)/1.0, 1)
           outfile = str(i)+'.xyz'
           f = open(outfile,"w")
           print >> f, no_of_atoms, "angstroem"
           print >> f, geometry, strained_lattice(lattice_a, i), strained_lattice(lattice_b, i), lattice_c
           for j in range(0,no_of_atoms):
               xyz_printing(atomic_symbol[j], strained_coordinate(x[j], i), strained_coordinate(y[j], i), z[j], mpol)
           pbs_printing(xyz_file_name, pbs_file_names)
       f.close( )

     
    elif strain_axis == 'xy-P':
      lattice_a_strain = lattice_a - strain_value_a
      for i1 in range(0,no_of_strain_interval):
        lattice_b_strain = lattice_b - strain_value_b
        for i in range(0,no_of_strain_interval):
           outfile = str(round(lattice_a_strain,3))+'-'+str(round(lattice_b_strain,3))+'.xyz'
           f = open(outfile,"w")
           print >> f, no_of_atoms, "angstroem"
           print >> f, geometry, lattice_a_strain, lattice_b_strain, lattice_c
           for j in range(0,no_of_atoms):
               x_strain[j] = x[j] / lattice_a * lattice_a_strain
               y_strain[j] = y[j] / lattice_b * lattice_b_strain
               xyz_printing(atomic_symbol[j], x_strain[j],  y_strain[j],  z[j], mpol)
           pbs_printing(xyz_file_name, pbs_file_names)
           lattice_b_strain = lattice_b_strain + (strain_value_b*2/no_of_strain_interval)
        lattice_a_strain = lattice_a_strain + (strain_value_a*2/no_of_strain_interval)
        f.close( ) 

    elif strain_axis == 'yz':
       for i in range(-4, 5):
           outfile = str(i)+'.xyz'
           f = open(outfile,"w")
           print >> f, no_of_atoms, "angstroem"
           print >> f, geometry, lattice_a, strained_lattice(lattice_b, i), strained_lattice(lattice_c, i)
           for j in range(0,no_of_atoms):
               xyz_printing(atomic_symbol[j], x[j], strained_coordinate(y[j], i), strained_coordinate(z[j], i), mpol)
           pbs_printing(xyz_file_name, pbs_file_names)
       f.close( )
     
    elif strain_axis == 'zx':
       for i in range(-4, 5):
           outfile = str(i)+'.xyz'
           f = open(outfile,"w")
           print >> f, no_of_atoms, "angstroem"
           print >> f, geometry, strained_lattice(lattice_a, i), lattice_b, strained_lattice(lattice_c, i)
           for j in range(0,no_of_atoms):
               xyz_printing(atomic_symbol[j], strained_coordinate(x[j], i), y[j], strained_coordinate(z[j], i), mpol)
           pbs_printing(xyz_file_name, pbs_file_names)
       f.close( )

     
    elif strain_axis == 'xyz-E':
       for i in range(-4, 5):
           i = round(float(i)/1.0, 1)
           outfile = str(i)+'.xyz'
           f = open(outfile,"w")
           print >> f, no_of_atoms, "angstroem"
           print >> f, geometry, strained_lattice(lattice_a, i), strained_lattice(lattice_b, i), strained_lattice(lattice_c, i)
           for j in range(0,no_of_atoms):
               xyz_printing(atomic_symbol[j], strained_coordinate(x[j], i), strained_coordinate(y[j], i), strained_coordinate(z[j], i), mpol)
           pbs_printing(xyz_file_name, pbs_file_names)
       f.close( )


     
    elif strain_axis == 'xyz-P':
      lattice_a_strain = lattice_a - strain_value_a
      for i1 in range(0,no_of_strain_interval):
        lattice_b_strain = lattice_b - strain_value_b
        for i in range(0,no_of_strain_interval):
            lattice_c_strain = lattice_c - strain_value_c
            for i2 in range(0,no_of_strain_interval):
                 outfile = str(round(lattice_a_strain,3))+'-'+str(round(lattice_b_strain,3))+'-'+str(round(lattice_c_strain,3))+'.xyz'
                 f = open(outfile,"w")
                 print >> f, no_of_atoms, "angstroem"
                 print >> f, geometry, lattice_a_strain, lattice_b_strain, lattice_c_strain
                 for j in range(0,no_of_atoms):
                     x_strain[j] = x[j] / lattice_a * lattice_a_strain
                     y_strain[j] = y[j] / lattice_b * lattice_b_strain
                     z_strain[j] = z[j] / lattice_c * lattice_c_strain
                     xyz_printing(atomic_symbol[j], x_strain[j],  y_strain[j],  z_strain[j], mpol)
                 pbs_printing(xyz_file_name, pbs_file_names)
                 lattice_c_strain = lattice_c_strain + (strain_value_c*2/no_of_strain_interval)
            lattice_b_strain = lattice_b_strain + (strain_value_b*2/no_of_strain_interval)
        lattice_a_strain = lattice_a_strain + (strain_value_a*2/no_of_strain_interval)
        f.close( )
    else:
      print "Error in typing proper input in strain axis"
      exit()
else:
    print "Error in typing apply strain(Y/N)"     
    exit()
############################################################################
                     # sub file  ###
############################################################################
sub_file_name = pbs_file_names
for i in range(0, size(sub_file_name)):
    dummy = sub_file_name[i]
    sub_file_name[i] = dummy[:-4]
f = open('sub.sh', 'w')
print >> f, "for i in ", ' '.join(sub_file_name)
print >> f, "do"
print >> f, "#qsub ./$i.pbs"
print >> f, "#rm -rf ./$i.pbs ./$i.xyz ./strain_inputs.txt"
print >> f, "#/W/dc200827/bigdft-suite/tmp/install/bin/bigdft-tool -n 64 --name=$i >$i.log"
print >> f, "done"
print "Total number of xyz-files created are : ", size(sub_file_name)
f.close()
f = open('strain_inputs.txt','w')
print >> f, "Strain applied on", strain_axis
print >> f, "The percentage of strain applied is -4 to +4"
print >> f, "No. of intervals are 9"
f.close()
cmd1 = "chmod +x sub.sh"
os.system(cmd1)
exit()
