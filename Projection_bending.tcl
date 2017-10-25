#This script is supposed to provide a formal method to calculate the principal
#axis of a selected part, especially a helix or CA of the hole protein.
#It will calculate the I matrix and its eigenvalue and eigenvector, choosing
#the smallest eigenvalue and its corresponding eigenvector as the pricipal 
#axis. Further it is going to calculate the polar and azimuthal angle



#Further modified to calculate the bendeing angle of projection
#on polar plane and azimuthal plane
#PROPERTY OF WENJUN ZHENG'S GROUP
#3/7/2015
#HAN WEN

#For simplicity and being realistic at the same time. I am going to only include#CA atoms, therefore please specify CA in selection, you can also modify for
#further use.

#Package la will be used, check for information on www.hume.com/la/
#Put la file in the same directory

package require La

namespace import -force La::*
set selection [atomselect top "resid 656 to 670 and chain D and name CA"]
set selection1 [atomselect top "resid 674 to 686 and chain D and name CA"]
set total_sel [atomselect top "resid 656 to 686 and chain D and name CA"]
#set selection [atomselect top "all and name CA"]
#modify for your choice
set mol_id [$selection molindex]
set angle_out [open "q_D_670-4.csv" w]
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
 

    #Already obtained principal axis, now do projection
    $total_sel frame $i
    set total_com [measure center $total_sel weight mass]
    #find the vector of com of total helix to define the plane
    #The normal vector for polar plane will be (-y,x,0)
    #The normal vector for azimuthal plane will be (x,y,0) 
    set cx [lindex $total_com 0]
    set cy [lindex $total_com 1]
    #Now calculate the projection on both plane:
    #V_projec =n*(V dot* n)/(n dot* n)     
    set para1_p [expr {($x1_v*(-$cy) + $y1_v*($cx))/($cx*$cx + $cy*$cy)}]
    set para2_p [expr {($x2_v*(-$cy) + $y2_v*($cx))/($cx*$cx + $cy*$cy)}]
    set para1_a [expr {($x1_v*($cx) + $y1_v*($cy))/($cx*$cx + $cy*$cy)}]
    set para2_a [expr {($x2_v*($cx) + $y2_v*($cy))/($cx*$cx + $cy*$cy)}]
    
    set x1_pn [expr {$para1_p*(-$cy)}] 
    set x2_pn [expr {$para2_p*(-$cy)}]
    set y1_pn [expr {$para1_p*$cx}]
    set y2_pn [expr {$para2_p*$cx}]

    set x1_an [expr {$para1_a*$cx}]
    set x2_an [expr {$para2_a*$cx}]
    set y1_an [expr {$para1_a*$cy}]
    set y2_an [expr {$para2_a*$cy}]

    #the projection on two planes
    set x1_p [expr {$x1_v - $x1_pn}]
    set y1_p [expr {$y1_v - $y1_pn}]
    set z1_p $z1_v
 
    set x2_p [expr {$x2_v - $x2_pn}]
    set y2_p [expr {$y2_v - $y2_pn}]
    set z2_p $z2_v

    set x1_a [expr {$x1_v - $x1_an}]
    set y1_a [expr {$y1_v - $y1_an}]
    set z1_a $z1_v

    set x2_a [expr {$x2_v - $x2_an}]
    set y2_a [expr {$y2_v - $y2_an}]
    set z2_a $z2_v


    #finally, calculate bending angle in both planes
    #normalize first
    set norm1p [expr {sqrt($x1_p * $x1_p + $y1_p*$y1_p + $z1_p * $z1_p)}]
    set x1_pnor [expr {$x1_p / $norm1p}]
    set y1_pnor [expr {$y1_p / $norm1p}]
    set z1_pnor [expr {$z1_p / $norm1p}]
    set norm2p [expr {sqrt($x2_p * $x2_p + $y2_p*$y2_p + $z2_p * $z2_p)}]
    set x2_pnor [expr {$x2_p / $norm2p}]
    set y2_pnor [expr {$y2_p / $norm2p}]
    set z2_pnor [expr {$z2_p / $norm2p}]
    
    set norm1a [expr {sqrt($x1_a * $x1_a + $y1_a*$y1_a + $z1_a * $z1_a)}]
    set x1_anor [expr {$x1_a / $norm1a}]
    set y1_anor [expr {$y1_a / $norm1a}]
    set z1_anor [expr {$z1_a / $norm1a}]
    set norm2a [expr {sqrt($x2_a * $x2_a + $y2_a*$y2_a + $z2_a * $z2_a)}]
    set x2_anor [expr {$x2_a / $norm2a}]
    set y2_anor [expr {$y2_a / $norm2a}]
    set z2_anor [expr {$z2_a / $norm2a}]

    set cos_theta_p [expr {$x1_pnor*$x2_pnor + $y1_pnor*$y2_pnor + $z1_pnor*$z2_pnor}]
    set theta_p [expr {180 * acos($cos_theta_p) / 3.1415926}]
   
    set cos_theta_a [expr {$x1_anor*$x2_anor + $y1_anor*$y2_anor + $z1_anor*$z2_anor}]
    set theta_a [expr {180 * acos($cos_theta_a) / 3.1415926}]


    puts $angle_out "$i,$theta_p,$theta_a"
    
}
close $angle_out
