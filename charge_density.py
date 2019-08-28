from numpy import *
from matplotlib.pyplot import *
from collections import *
from operator import itemgetter
#############################################################################
                  # Reading charge density
#############################################################################
def Read_charge(file_name):
   print "Reading", file_name, "..."
   charge = []
   i = 0
   with open(file_name, 'r') as position_data:
      for position_line in position_data:
          words = position_line.split()
          if i == 0:
             dummy = 0
          elif i == 1:
             dummy = 0
          elif i == 2:
             no_of_atoms = int(words[0])
          elif i == 3:
             Nx, Xvoxel = int(words[0]), float(words[1])
          elif i == 4:
             Ny, Yvoxel = int(words[0]), float(words[2])
          elif i == 5:
             Nz, Zvoxel = int(words[0]), float(words[3])
          elif i > no_of_atoms + 5:   
               for j in range(0, size(words)):
                   charge.append(float(words[j])*Xvoxel*Yvoxel*Zvoxel)
          i += 1
   tot_charge = sum(charge)
   lattice_a   = Nx*Xvoxel
   lattice_b   = Ny*Yvoxel
   lattice_c   = Nz*Zvoxel 
   return Xvoxel, Yvoxel, Zvoxel, lattice_a, lattice_b, lattice_c, tot_charge, charge, Nx, Ny, Nz
################################################################################
def radius(x_f, y_f, z_f):
    return(float(round(sqrt((x_f)**2 + (y_f)**2 + (z_f)**2),2)))

################################################################################
                        # FINDING CHARGE DENSITY
################################################################################
def charge_desity(initial_x, initial_y, initial_z, charge_f,  Xvoxel, Yvoxel, Zvoxe, Nx, Ny, Nz):
   II1 = int(initial_x/0.529177/Xvoxel)
   II2 = int(initial_y/0.529177/Yvoxel)
   II3 = int(initial_z/0.529177/Zvoxel) 
   r = []
   for i in range(0, Nx):
       d1 = (i-II1)*Xvoxel*0.529177
       for j in range(0, Nx):
           d2 = (j-II2)*Xvoxel*0.529177
           for k in range(0, Nx):
               d3 = (k-II3)*Xvoxel*0.529177
               r.append(radius(d1, d2, d3))
   
   r_and_charge = transpose(sorted(transpose([r, charge_f]), key=itemgetter(0)))
   charge = []
   r = []
   r = r_and_charge[0]
   charge = r_and_charge[1]
   
   r_sort_f = []
   charge_sort_f = [] #consicutive addition of charge
   
   i = 0
   dummy1 = 0
   dummy2 = 0
   while(i < size(r)-1):
       d1 = charge[i]
       while(r[i+1] == r[i]):
           d1 += charge[i]
           i += 1
       else:
           r_sort_f.append(r[i])
           dummy1 += d1
           charge_sort_f.append(dummy1)
       i += 1
   if r[size(r)-1] != r[size(r)-2]:
       r_sort_f.append(r[size(r)-1])
       x1 = dummy1+charge[size(r)-1]
       charge_sort_f.append(x1)
   return r_sort_f, charge_sort_f

################################################################################
output                          = Read_charge("electronic_density.cube")
Xvoxel, Yvoxel, Zvoxel          = output[0], output[1], output[2]
lattice_a, lattice_b, lattice_c = output[3], output[4], output[5]
tot_charge, charge              = output[6], output[7]
Nx, Ny, Nz                      = output[8], output[9], output[10]

print
print "Total charges: ", tot_charge
print

i_c = [] # input_coordinates = [], contains name of xyz coordinate, x, y, z
no_of_centers = 0 # centers to calculate charge densities
with open('input_coordinates.txt', 'r') as position_data:
      for position_line in position_data:
          no_of_centers +=1
          w = position_line.split()
          w[0]= str(w[0])
          w[1]= float(w[1])
          w[2]= float(w[2])
          w[3]= float(w[3])
          i_c.append([w[0], w[1], w[2], w[3]])

x =[]
y =[]
for j in range(0, no_of_centers):
    cmd = "Calculating charge density at "+i_c[j][0]+"..."
    print cmd
    print
    r_sort = []
    charge_sort = []
    r_sort, charge_sort = charge_desity(i_c[j][1], i_c[j][2], i_c[j][3], charge, Xvoxel, Yvoxel, Zvoxel, Nx, Ny, Nz)
    x_value = []
    y_value = []
    for i in range(0, size(r_sort)):
       if r_sort[i] > 0.5 and r_sort[i] < 2.5:
          x_value.append(r_sort[i])
          y_value.append(charge_sort[i])
    x.append(x_value)
    y.append(y_value)
    #y_value = diff(y_value)
    #del(x_value[0])

figure(figsize = (5, 8))
for i in range(0,3):
    labels = "From "+i_c[i][0]
    plot(x[i], y[i], label=labels)
    xlabel("Radial distance from the vacancy ($\AA$)")
    ylabel("Number of electrons")
    legend(loc = 'upper left')#
    grid(True)
savefig('radial_charges_1', bbox_inches='tight') 

figure(figsize = (5, 8))
for i in range(3,6):
    labels = "From "+i_c[i][0]
    plot(x[i], y[i], label=labels)
    xlabel("Radial distance from the vacancy ($\AA$)")
    ylabel("Number of electrons")
    legend(loc = 'upper left')#
    grid(True)
savefig('radial_charges_2', bbox_inches='tight')

figure(figsize = (5, 8))
for i in range(6,9):
    labels = "From "+i_c[i][0]
    plot(x[i], y[i], label=labels)
    xlabel("Radial distance from the vacancy ($\AA$)")
    ylabel("Number of electrons")
    legend(loc = 'upper left')#
    grid(True)
savefig('radial_charges_3', bbox_inches='tight')
figure(figsize = (5, 8))
for i in range(9,12):
    labels = "From "+i_c[i][0]
    plot(x[i], y[i], label=labels)
    xlabel("Radial distance from the vacancy ($\AA$)")
    ylabel("Number of electrons")
    legend(loc = 'upper left')#
    grid(True)
savefig('radial_charges_4', bbox_inches='tight')

print
print "Done"
