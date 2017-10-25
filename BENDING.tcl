#This script is supposed to provide a formal method to calculate the principal
#axis of a selected part, especially a helix or CA of the hole protein.
#It will calculate the I matrix and its eigenvalue and eigenvector, choosing
#the smallest eigenvalue and its corresponding eigenvector as the pricipal 
#axis. Further it is going to calculate the polar and azimuthal angle
#PROPERTY OF WENJUN ZHENG'S GROUP
#2/10/2015
#HAN WEN

#For simplicity and being realistic at the same time. I am going to only include#CA atoms, therefore please specify CA in selection, you can also modify for
#further use.

#Package la will be used, check for information on www.hume.com/la/
#Put la file in the same directory

package require La

namespace import -force La::*
set selection [atomselect top "resid 656 to 665 and chain D and name CA"]
set selection1 [atomselect top "resid 666 to 686 and chain D and name CA"]
#set selection [atomselect top "all and name CA"]
#modify for your choice
set mol_id [$selection molindex]
set angle_out [open "s6wjt330bending_D.csv" w]
set num_frames [molinfo $mol_id get numframes]
puts $num_frames
for {set i 0} {$i < $num_frames} {incr i} {
#Loop for all frames
    $selection frame $i
    set x_tem [$selection get x]
    set y_tem [$selection get y]
    set z_tem [$selection get z]
    #get a list of all x, y ,z
    set I_xx 0
    set I_xy 0
    set I_xz 0
    set I_yy 0
    set I_yz 0
    set I_zz 0
    set com [measure center $selection weight mass]
    #calculate the center of mass
    #The components of inertia matrix
    foreach xi $x_tem yi $y_tem zi $z_tem {
	set xi [expr {$xi - [lindex $com 0]}]
	set yi [expr {$yi - [lindex $com 1]}]
	set zi [expr {$zi - [lindex $com 2]}]
 	set I_xx [expr {$I_xx + $yi*$yi + $zi*$zi}]
        set I_xy [expr {$I_xy - $xi*$yi}]
        set I_xz [expr {$I_xz - $xi*$zi}]
        set I_yy [expr {$I_yy + $xi*$xi + $zi*$zi}]
	set I_yz [expr {$I_yz - $yi*$zi}]
	set I_zz [expr {$I_zz + $xi*$xi + $yi*$yi}]
    }
    set I_matrix [list 2 3 3 $I_xx $I_xy $I_xz $I_xy $I_yy $I_yz $I_xz $I_yz $I_zz]
    #puts [show $I_matrix]
    #Set the inertia matrix, 2 means 2-D ,3 Row, 3 Column, therefore lindex strart at 3
    mevsvd_br I_matrix eigenvalue
    #La packet function, return I with column as eigenvectors eigenvalue in decreasing order
    set principal_axis "[lindex $I_matrix 5] [lindex $I_matrix 8] [lindex $I_matrix 11]"
    #puts $principal_axis
    #puts [show $eigenvalue]
    set x1_v [lindex $I_matrix 5]
    set y1_v [lindex $I_matrix 8]
    set z1_v [lindex $I_matrix 11]
    #Calculate two angles, same as simple version, polar one is between it and z
    #azimuthal one is the projection on xy plane with x








    $selection1 frame $i   
    set x_tem1 [$selection1 get x]
    set y_tem1 [$selection1 get y]
    set z_tem1 [$selection1 get z]
    #get a list of all x, y ,z
    set I_xx 0
    set I_xy 0
    set I_xz 0
    set I_yy 0
    set I_yz 0
    set I_zz 0
    set com1 [measure center $selection1 weight mass]
    #calculate the center of mass
    #The components of inertia matrix
    foreach xi $x_tem1 yi $y_tem1 zi $z_tem1 {
        set xi [expr {$xi - [lindex $com1 0]}]
        set yi [expr {$yi - [lindex $com1 1]}]
        set zi [expr {$zi - [lindex $com1 2]}]
        set I_xx [expr {$I_xx + $yi*$yi + $zi*$zi}]
        set I_xy [expr {$I_xy - $xi*$yi}]
        set I_xz [expr {$I_xz - $xi*$zi}]
        set I_yy [expr {$I_yy + $xi*$xi + $zi*$zi}]
        set I_yz [expr {$I_yz - $yi*$zi}]
        set I_zz [expr {$I_zz + $xi*$xi + $yi*$yi}]
    }
    set I_matrix [list 2 3 3 $I_xx $I_xy $I_xz $I_xy $I_yy $I_yz $I_xz $I_yz $I_zz]
    #puts [show $I_matrix]
    #Set the inertia matrix, 2 means 2-D ,3 Row, 3 Column, therefore lindex strart at 3
    mevsvd_br I_matrix eigenvalue
    #La packet function, return I with column as eigenvectors eigenvalue in decreasing order
    set principal_axis "[lindex $I_matrix 5] [lindex $I_matrix 8] [lindex $I_matrix 11]"
    #puts $principal_axis
    #puts [show $eigenvalue]
    set x2_v [lindex $I_matrix 5]
    set y2_v [lindex $I_matrix 8]
    set z2_v [lindex $I_matrix 11]
 

    
    set norm1 [expr {sqrt($x1_v * $x1_v + $y1_v*$y1_v + $z1_v * $z1_v)}]
    set x1_n [expr {$x1_v / $norm1}]
    set y1_n [expr {$y1_v / $norm1}]
    set z1_n [expr {$z1_v / $norm1}]
    set norm2 [expr {sqrt($x2_v * $x2_v + $y2_v*$y2_v + $z2_v * $z2_v)}]
    set x2_n [expr {$x2_v / $norm2}]
    set y2_n [expr {$y2_v / $norm2}]
    set z2_n [expr {$z2_v / $norm2}]
    
    
    set cos_theta [expr {$x1_n*$x2_n + $y1_n*$y2_n + $z1_n*$z2_n}]
    set theta [expr {180 * acos($cos_theta) / 3.1415926}]
    puts $angle_out "$theta,$i"
    
}
close $angle_out
