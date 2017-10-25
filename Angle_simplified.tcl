#This program is supposed to measure the angles of a helix over multi-frames
#Using a simplified model, assume Z axis is the symmetrical axis
#the helix using 3D least square fitting
#PROPERTY OF WENJUN ZHENG'S GROUP
#2/5/2015
#HAN WEN
#This script can be used for analysis of helix as a simple 'peek', for it is sequence depending.



#Assuming you already load the trajectory (of lipid protein system) in pdb format
#This script is supposed to deal with gromacs trj, therefore used chain A in selection
set selection [atomselect top "resid 416 to 428 and chain A and name CA"]
#set selection [atomselect top "all and name CA and z > 170"]
#Modify the selection to deal with different parts

#Procedure to calculate least square fit in each dimension, by fitting to x = x_out*i + b
# x_out = 12/( (N(N^2 - 1)) ) sum[ (i-(N-1)/2) * xi]
#proc lsq { x } {
#    set num_list [llength $x]
#    set x_out 0                        
#    set n 0.0                  
#    set d [expr {0.5 * ($num_list - 1)}]
#    foreach coor $x {
#	set x_out [expr {$x_out + ($n - $d) * $coor}]
#	set n [expr {$n + 1.0}]
#    }
#    return $x_out
#}
proc lsq { x } {
  set N [llength $x]
  set xtot 0
  #puts $x
  set d [expr {0.5*($N - 1)}]
  
  set i 0.0
  foreach elem $x {
    set xtot [expr {$xtot + ($i - $d) * $elem }]
    set i [expr {$i + 1.0}]
  }
  
# no need to normalize if all we want is the direction
  set xtot [expr $xtot * 12 / ($N * ($N * $N - 1))]
  return $xtot
}    



set mol_id [$selection molindex]
#Get selected index for further use

set angle_out [open "angletry" w]

set num_frames [molinfo $mol_id get numframes] 
#Get total frames number

for {set i 0} {$i < $num_frames} {incr i} {
#The loop for all frames
    $selection frame $i
#    lappend com [measure center $selection weight mass]
    #calculate the center of mass
    set x_tem [$selection get x]        
    set y_tem [$selection get y]
    set z_tem [$selection get z]
    # x_v the x component for result fitting vector
    #puts $x_tem
    #puts $y_tem
    #puts $z_tem
    set x_v [lsq $x_tem] 
    set y_v [lsq $y_tem]
    set z_v [lsq $z_tem]
    #puts $x_v
    #Normalization
    set norm [expr {sqrt($x_v * $x_v + $y_v * $y_v +  $z_v * $z_v)}]
    set x_n [expr {$x_v / $norm}]
    set z_n [expr {$z_v / $norm}]
    set y_n [expr {$y_v / $norm}]
    #Calculate theta_p and theta_a
    puts "vec $x_v; $y_v ; $z_v"
    set norm_p [expr {sqrt($x_v * $x_v + $z_v * $z_v)}]   
    set z_p [expr {$z_v / $norm_p}]
    set cos_theta_p $z_p
    set norm_a [expr {sqrt($x_v * $x_v + $y_v * $y_v)}]
    set x_a [expr {$x_v / $norm_a}]
    set cos_theta_a $x_a
    set theta_p [expr {180 * acos($cos_theta_p) / 3.1415926}]
    set theta_a [expr {180 * acos($cos_theta_a) / 3.1415926}]
    puts $angle_out "$theta_p ; $theta_a"
} 

#puts $com   
