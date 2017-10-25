#This script is supposed to provide a formal method to calculate the principal
#axis of a selected part, especially a helix or CA of the protein.
#It will calculate the I matrix and its eigenvalue and eigenvector, choosing
#the smallest eigenvalue and its corresponding eigenvector as the pricipal 
#axis. Further it is going to calculate the polar and azimuthal angle
#by project into plane formed by z axis and COM
#PROPERTY OF WENJUN ZHENG'S GROUP
#5/3/2015
#HAN WEN

#For simplicity and being realistic at the same time. I am going to only include#CA atoms, therefore please specify CA in selection, you can also modify for
#further use.

#Package la will be used, check for information on www.hume.com/la/
#Put la file in the same directory

package require La

namespace import -force La::*
set selection [atomselect top "resid 630 to 642 and segname PROA and name CA"]
#modify for your choice
set mol_id [$selection molindex]
                       set angle_out [open "529tilt_pore_A.csv" w]
set num_frames [molinfo $mol_id get numframes]
puts $num_frames

$selection frame 0
set total_com [measure center $selection weight mass]
set cx_o [lindex $total_com 0]
set cy_o [lindex $total_com 1]
set cz_o [lindex $total_com 2]
set zdiff [expr {0}]
set theta_com 0
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
    
    #now do projection
    #set total_com [measure center $selection weight mass]

    set cx [lindex $com 0]
    set cy [lindex $com 1]
    set cz [lindex $com 2]

    #The normal vector for polar plane will be (-y,x,0)
    #The normal vector for azimuthal plane will be (x,y,0)
    #Now calculate the projection on both plane:
    #V_projec =n*(V dot* n)/(n dot* n)
    set para_p [expr {($x_v*(-$cy) + $y_v*($cx))/($cx*$cx + $cy*$cy)}]
    set para_a [expr {($x_v*($cx) + $y_v*($cy))/($cx*$cx + $cy*$cy)}]

    set x_pn [expr {$para_p*(-$cy)}]
    set y_pn [expr {$para_p*$cx}]

    set x_an [expr {$para_a*$cx}]
    set y_an [expr {$para_a*$cy}]

    set x_p [expr {$x_v - $x_pn}]
    set y_p [expr {$y_v - $y_pn}]
    set z_p $z_v

    set x_a [expr {$x_v - $x_an}]
    set y_a [expr {$y_v - $y_an}]
    set z_a $z_v

    #Now do normalization
    set normp [expr {sqrt($x_p * $x_p + $y_p*$y_p + $z_p * $z_p)}]
    #set x_pnor [expr {$x_p / $normp}]
    #set y_pnor [expr {$y_p / $normp}]
    set z_pnor [expr {$z_p / $normp}]

    set norma [expr {sqrt($x_a * $x_a + $y_a*$y_a + $z_a * $z_a)}]
    #set x_anor [expr {$x_a / $norma}]
    #set y_anor [expr {$y_a / $norma}]
    set z_anor [expr {$z_a / $norma}]

    set cos_theta_p $z_pnor
    set cos_theta_a $z_anor

    set theta_p [expr {180 * acos($cos_theta_p) / 3.1415926}]
    set theta_a [expr {180 * acos($cos_theta_a) / 3.1415926}]
    
    #Calculate movement of COM: 
    #azimuthal angle:
    #puts "$cx,$cy,$cz,$cx_o,$cy_o"
    if {$i == 0} {
        set xpref [expr {$x_p}]       
        set xaref [expr {$x_a}]  
    }

    if {$i>0} {
        set cos_com [expr {($cx*$cx_o + $cy*$cy_o)/ ( (sqrt($cx*$cx+$cy*$cy))*(sqrt($cx_o*$cx_o+$cy_o*$cy_o)))}]
        set theta_com [expr {180 * acos($cos_com) / 3.1415926}]
        set zdiff [expr {$cz-$cz_o}]

        if { $cx*$cy_o-$cy*$cx_o <0 } {
            set theta_com [expr {-$theta_com}]      
        } 
#	if { $xpref*$x_p <=0 } {
#	    set theta_p [expr {-$theta_p}]
#	}
#    	if { $xaref*$x_a <=0 } {
#            set theta_a [expr {-$theta_a}]
#	}
    }
    puts $angle_out "$i,$theta_p , $theta_a,$theta_com,$zdiff, xyz:,  $cx,$cy,$cz"
 
}
close $angle_out
