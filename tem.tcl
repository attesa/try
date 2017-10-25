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
set selection [atomselect top "resid 630 to 642 and chain D and name CA"]
#set selection [atomselect top "all and name CA"]
#modify for your choice
set mol_id [$selection molindex]
set angle_out [open "all_pore_align_chain_D.dat" w]
set num_frames [molinfo $mol_id get numframes]

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
    set x_v [lindex $I_matrix 5]
    set y_v [lindex $I_matrix 8]
    set z_v [lindex $I_matrix 11]
    #Calculate two angles, same as simple version, polar one is between it and z
    #azimuthal one is the projection on xy plane with x
    if "$y_v < 0.0" {
	set x_v [expr {-$x_v}]
	set y_v [expr {-$y_v}]
	set z_v [expr {-$z_v}]
    } 
    #Alternate the direction to ensure theta_a from 0 to 180
    #comment four lines for pore



    
    set norm [expr {sqrt($x_v * $x_v + $y_v*$y_v + $z_v * $z_v)}]
    set z_n [expr {$z_v / $norm}]
    set cos_theta_p $z_n
    set theta_p [expr {180 * acos($cos_theta_p) / 3.1415926}]

    set norm_a [expr {sqrt($x_v * $x_v + $y_v * $y_v)}]
    set x_a [expr {$x_v / $norm_a}]
    set cos_theta_a $x_a
    #puts $principal_axis
    #if "$y_v < 0.0" {
#	set theta_a [expr {360-180 * acos($cos_theta_a) / 3.1415926}]
#    } else {
    set theta_a [expr {180 * acos($cos_theta_a) / 3.1415926}]
#    }
#uncomments for pore
    puts $angle_out "$theta_p ; $theta_a"
 
}
close $angle_out
