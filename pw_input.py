from numpy import *
from matplotlib.pyplot import *
from collections import *
from operator import itemgetter
import sys, string, os, commands
##########################################################################

##########################################################################
print "******************************************"
print "            The input options:"
print "******************************************"
print "1. Use '.vasp' file or"
print "2. Use '.out' QE file"
input_option = input(" Type 1 or 2: ", )
print "******************************************"
##########################################################
            #  Reading .vasp file
###########################################################
if input_option == 1:
    vasp_data = []
    input_files=raw_input("Write the name of input file(without .vasp):  ")
    vasp_file_name = input_files+'.vasp'
    #geometry = raw_input("Is it free, surface, or periodic:  ")
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
    lattice_a = lattice_a/0.529177249
    lattice_b = lattice_b/0.529177249
    lattice_c = lattice_c/0.529177249
   
    types_of_atoms = vasp_data[5]
    no_of_types_of_atoms = size(types_of_atoms)
    
    no_of_atoms_in_each_type = vasp_data[6]
    for i in range(0, no_of_types_of_atoms):
        no_of_atoms_in_each_type[i] = int(no_of_atoms_in_each_type[i])  # to convert int
    no_of_atoms = sum(no_of_atoms_in_each_type)
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
    scale = ''.join(vasp_data[7])
    if scale in ('Direct', 'direct'):
       for i in range(0, no_of_atoms):
            x.append(xyz_data[i,0])
            y.append(xyz_data[i,1])
            z.append(xyz_data[i,2])
    else:
        for i in range(0, no_of_atoms):
            x.append(xyz_data[i,0]/0.529177249)
            y.append(xyz_data[i,1]/0.529177249)
            z.append(xyz_data[i,2]/0.529177249)
    atomic_symbol = []
    for i in range(0, no_of_types_of_atoms):
          for j in range(0,no_of_atoms_in_each_type[i]):
                atomic_symbol.append(types_of_atoms[i])
    

#for i in range(0,sum(no_of_atoms_in_each_type)):
#    if y[i] > lattice_b-1:
#          y[i] = y[i] - lattice_b
#    else:
#          y[i] = y[i]

    print "ibrav values:"
    print "0 :free (crystal axis provided in input card CELL_PARAMETERS)"
    print "1: cubic P (sc), v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,1)"
    print "2: cubic F (fcc), v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)"
    print "3: cubic I (bcc), v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)"
    print "4: Hexagonal and Trigonal P, celldm(3)=c/a, v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)"
    print "8: Orthorhombic P, celldm(2)=b/a, celldm(3)=c/a, v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)"
    ibrav = input("ibvrav value:")

############################################################
          # Reading .out QE file
############################################################
elif input_option == 2:
     input_files = raw_input("Write name of 'output' file without '.out': ")
     #geometry = raw_input("Is it 'surface' or 'periodic' structure?: ")
     qefile = input_files+".out"
     cmd1 = "grep -e 'number of atomic types' "+qefile+" | awk '{print $6}' | tail -n1"
     status, output = commands.getstatusoutput(cmd1)
     no_of_types_of_atoms = int(output)
     cmd2 = "grep -e 'number of atoms/cell' "+qefile+" | awk '{print $5}' | tail -n1"
     status, output = commands.getstatusoutput(cmd2)
     no_of_atoms = int(output)
     cmd3 = "grep -e ' bravais-lattice index' "+qefile+" | awk '{print $4}' | tail -n1"
     status, output = commands.getstatusoutput(cmd3)
     ibrav = int(output)
     cmd4 = "grep -e 'celldm(1)' "+qefile+" | tail -n1 | awk '{print $2}' | tail -n1"
     status, output = commands.getstatusoutput(cmd4)
     celldm1 = float(output)
     cmd5 = "grep -e 'celldm(1)' "+qefile+" | tail -n1 | awk '{print $4}' | tail -n1"
     status, output = commands.getstatusoutput(cmd5)
     celldm2 = float(output)
     cmd6 = "grep -e 'celldm(1)' "+qefile+" | tail -n1 | awk '{print $6}' | tail -n1"
     status, output = commands.getstatusoutput(cmd6)
     celldm3 = float(output)

     lattice_a   = float(celldm1)
     lattice_b   = float(celldm2)*lattice_a
     lattice_c   = float(celldm3)*lattice_a

     cmd7 = "grep -A"+str(no_of_atoms)+" 'ATOMIC_POSITIONS' "+qefile+" | tail -n"+str(no_of_atoms)
     status, output = commands.getstatusoutput(cmd7)
     position_data = output
    
     atomic_symbol = []
     x = []
     y = []
     z = []
    
     position_line = position_data.splitlines()
     for i in range(0, no_of_atoms):
         dummy = position_line[i]
         dummy = dummy.split()
         atomic_symbol.append(str(dummy[0]))

         cmd8 = "grep -e 'ATOMIC_POSITIONS' "+qefile+" | tail -n1"  
         status, output = commands.getstatusoutput(cmd8)
         scale = 'angstrom'
         if scale in output:
             x.append(float(dummy[1]))
             y.append(float(dummy[2]))
             z.append(float(dummy[3]))
         else:        
             x.append(float(dummy[1])*lattice_a)
             y.append(float(dummy[2])*lattice_b)
             z.append(float(dummy[3])*lattice_c)

     types_of_atoms = list(set(atomic_symbol))
     no_of_atoms_in_each_type = []
     for i in range (0, no_of_types_of_atoms):
         j = 0
         for k in range(0, no_of_atoms):
             if types_of_atoms[i] == atomic_symbol[k]:
                j += 1
         no_of_atoms_in_each_type.append(str(j))
############################################################################
                 #Pseudopoetentials
############################################################################
pseudopotentials = ['Cd 112.411  Cd.pz-dn-kjpaw_psl.0.2.UPF', 'Te 127.600  Te.pz-dn-kjpaw_psl.0.2.2.UPF']
def potential(fx):
    if fx in ''.join(pseudopotentials):
        for i in range(size(pseudopotentials)):
              if fx in pseudopotentials[i]:
                 fy = pseudopotentials[i]
                 return fy
    else:
        print fx+"-pseudopotential is not available"
        exit()
print "Available pseudopotentials:"
for i in range(0, no_of_types_of_atoms):
    print potential(types_of_atoms[i])
print "*******************************"
############################################################################
calculation = raw_input("vc-relax, relax, scf, nscf, or bands:")
calculation = "\'"+calculation+"\'"
#print "###################################"

#########################################################################

f = open(input_files+".in", "w")
#############################
      ## &CONTROL ###
#############################
print >> f, "&CONTROL"
print >> f, "   calculation    =", calculation  
print >> f, "   verbosity      = 'high'"
print >> f, "   restart_mode   = 'from_scratch'"
print >> f, "   wf_collect     = .FALSH."
print >> f, "   nstep          =  200"
#print >> f, "   tstress        = ", tstress
#print >> f, "   tprnfor        = ", tprnfor
print >> f, "   outdir         = '/W/ss256439/QE/temp'" 
print >> f, "   pseudo_dir     = '/home/ss256439/myopt/QE-pot/'"
print >> f, "   prefix         =","\'"+''.join(types_of_atoms)+"-"+input_files+"\'"
print >> f, "   etot_conv_thr  =  1.0D-5"
print >> f, "   forc_conv_thr  =  1.0D-4"
               
print >> f, "/"
#############################
     ## &SYSTEM ##
#############################
nat = no_of_atoms
ntyp = no_of_types_of_atoms
spin_atoms = ['Mn', 'Fe', 'Co', 'Ni']
if spin_atoms in types_of_atoms:
   nspin = 2
else:
   nspin = 1

print >> f, "&SYSTEM"
print >> f, "   ibrav          = ", ibrav
if ibrav == 1:
   print >> f, "   celldm(1)      = ", lattice_a
elif ibrav == 4:
   print >> f, "   celldm(1)      = ", lattice_a
   print >> f, "   celldm(3)      = ", lattice_c/lattice_a
elif ibrav == 8:
   print >> f, "   celldm(1)      = ", lattice_a
   print >> f, "   celldm(2)      = ", lattice_b/lattice_a
   print >> f, "   celldm(3)      = ", lattice_c/lattice_a

print >> f, "   nat            = ", nat
print >> f, "   ntyp           = ", ntyp
#print >> f, "   nbnd           = ", nbnd
print >> f, "   nspin          = ", nspin
if nspin == 2:
   for i in range(0, ntyp):
       print f >> "   starting_magnetization("+str(i)+") =  0.3"

ecutwfc = 45   # input("ecut value (for 111--> 50) (for 001--> 40): ")
ecutrho = 400  #ecutwfc*6
nosym   = ".FALSE."
occupations = "'smearing'"
degauss = 0.001  #input("degauss for 111--> 0.005 anf for 001--> 0.01: ")
smearing  = "'gaussian'"
print >> f, "   ecutwfc        = ", ecutwfc
print >> f, "   ecutrho        = ", ecutrho
print >> f, "   nosym          =", nosym
print >> f, "   occupations    =", occupations
print >> f, "   degauss        = ", degauss
print >> f, "   smearing       =", smearing
print >> f, "/"
#############################
     ## &ELECTRONS ##
#############################
print >> f, "&ELECTRONS" 

print >> f, "   electron_maxstep   =  100"
print >> f, "   conv_thr           =  1.D-5"
mixing_mode = "'plain'"
mixing_beta = 0.7
diagonalization = "'david'"
print >> f, "   mixing_mode        =", mixing_mode
print >> f, "   mixing_beta        = ", mixing_beta
print >> f, "   diagonalization    =", diagonalization
print >> f, "/"
#############################
     ## &IONS ##
#############################
if calculation in ("'relax'", "'vc-relax'"):
    print >> f, "&IONS"
    dynamics = "'bfgs'"
    print >> f, "   ion_dynamics  =", dynamics
    print >> f, "/"
#############################
     ## &CELL ##
#############################
if calculation == "'vc-relax'":
    print >> f, "&CELL"
    print >> f, "   cell_dynamics =", dynamics
    print >> f, "   press         =  0.D0"
    print >> f, "   cell_factor   =  2.0"
    print >> f, "   cell_dofree   = 'all'"
    print >> f, "/"
#############################
     ## ATOMIC_SPECIES ##
#############################
print >> f, "ATOMIC_SPECIES"
for i in range(0, ntyp):
      print >> f, potential(types_of_atoms[i]) #, mass(i), pp_file_name(i)
#############################
     ## K_POINTS { tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c } ##
#############################
print >> f, "K_POINTS {automatic}"
if calculation in["'vc-relax'", "'relax'", "scf'"]:
    print "No. of k-points along x, y, z:"
    kx = input("kx = ",)
    ky = input("ky = ",)
    kz = input("kz = ",)
    if kx > 1:
       kx_shift = 1
    else:
       kx_shift = 0
    if ky > 1:
       ky_shift = 1
    else:
       ky_shift = 0
    if kz > 1:
       kz_shift = 1
    else:
       kz_shift = 0
    print >> f, kx, ky, kz, kx_shift, ky_shift, kz_shift

#############################
     ## ATOMIC_POSITIONS { alat | bohr | angstrom | crystal | crystal_sg } ##
#############################
#print >> f, "ATOMIC_POSITIONS { angstrom }"
if scale in ('Direct', 'direct'):
   print >> f, "ATOMIC_POSITIONS { crystal }"
   for k in range(0, no_of_atoms):
       print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(atomic_symbol[k], x[k],  y[k],  z[k]))
else:
   print >> f, "ATOMIC_POSITIONS { alat }"
   for k in range(0, no_of_atoms):
     x[k] = (x[k])/lattice_a 
     y[k] = (y[k])/lattice_a 
     z[k] = (z[k])/lattice_a
     print >> f, ("%-6s %11.6f    %11.6f    %11.6f" %(atomic_symbol[k], x[k],  y[k],  z[k]))

print "********************************"
print "      1. "+input_files+".in ", "is created!"

#############################################
    # PBS file
############################################

f = open(input_files+".pbs", "w")
print >> f, "#!/bin/bash -l"
print >> f, "#PBS -q csp"
print >> f, "#PBS -N "+input_files
print >> f, "#PBS -j oe"
print >> f, "#PBS -l nodes=2:ppn=32"
print >> f, "#PBS -l walltime=24:00:00"
print >> f, "setenv MODULEPATH \"/home/prog/usr/prod/modulefiles/tools/:/home/prog/usr/prod/modulefiles/applications/:/home/prog/usr/prod/modulefiles/libraries/\""
print >> f, "module spider espresso/6.3 mkl/18"
print >> f, "cd $PBS_O_WORKDIR"
print >> f, "limit stacksize unlimited"
print >> f, "mpirun /home/prog/usr/prod/applications/espresso/6.3/bin/pw.x -in "+input_files+".in > "+input_files+".out"
print "      2. "+input_files+".PBS file also created!"
print "********************************"
