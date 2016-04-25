#!/bin/bash

grid=1080


##########################
for i in $(echo 10xyz); do
count=10
##########################
#
#
##############################################################
#count=2
#for i in $(echo 3xyz 4xyz 5xyz 6xyz 7xyz 8xyz 9xyz 10xyz); do
#count=$((count+1))
##############################################################
#
#
cat > input.$count << EOF
inputhelix NewNievergelt.$i
helixout_name testhelix.$count.$grid.out
coord_type 0
num_grid $grid
natoms $count
nframes 1
grid_phi_beg 90
grid_phi_end 135
grid_theta_beg 90
grid_theta_end 135
print_to_plot 1
print_step 0
EOF
#
#
/Users/k/Desktop/Manuscripts/kHelix/Code/kHelix.o input.$count
#
#
done
