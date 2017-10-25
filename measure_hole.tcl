#This script is supposed to measure the hole in membrane protein.
#To be specific, it is going to measure the target resid and find the minimal
#hole radius through out the trajectory
#To be safe, it will use all atoms not just CA, yet this setting can be easily
#modified in atomselection part, shown below

#To use this script, the hole.tcl from vmd website is needed in the same
#directory and the software hole need to be installed, then change the variable
#in the hole.tcl so that it can be linked to hole software.

#Align about s1-s6 first with p before use

#PROPERTY OF WENJUN ZHENG'S GROUP
#2/18/2015 (Chinese new year)
#HAN WEN

source hole.tcl

set selection [atomselect top protein]
set output [open "holeout_t330_115res679.dat" w]
#set output file 
set mol_id [$selection molindex]
set num_frames [molinfo $mol_id get numframes]
for {set i 0} {$i < $num_frames} {incr i} {
#Loop for all frames
    $selection frame $i
    puts $i
    set hole_list [::Hole::runhole $selection]
    #hole_list contains all the returned value in form of 
    #[z radius resname resid]
    set target_resid 679
    #The resid you want to analysis
    set r_out 999
    #Out_put radius
    set z_out -1
    #Out_put z coordinate
    foreach element $hole_list {
	set resid_in_list [lindex $element 3]
        if {$resid_in_list == $target_resid} then {
	    set r_tem [lindex $element 1]
            if { $r_tem < $r_out} then {
		set r_out $r_tem
		set z_out [lindex $element 0]
	    }
	}
    }
    puts $output "$i;$r_out;$z_out"
}
close $output

