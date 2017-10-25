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

set selection [atomselect top "name BB"]
set output [open "holec5t330_mar_673.csv" w]
set output1 [open "holec5t330_mar_679.csv" w]
set output2 [open "holec5t330_mar_643.csv" w]
set outputg [open "holec5t330_mar_global.csv" w]
#set output file 
#set output [open "r671" w]
#set output1 [open "r679.dat" w]
#set output2 [open "r643.dat" w]
#set outputg [open "g.dat" w]



set mol_id [$selection molindex]
set num_frames [molinfo $mol_id get numframes]
for {set i 0} {$i < $num_frames} {incr i 4} {
#Loop for all frames
    $selection frame $i
    puts $i
    set hole_list [::Hole::runhole $selection]
    #hole_list contains all the returned value in form of 
    #[z radius resname resid]
    #set target_resid 673
    #set target_resid1 679
    #set target_resid2 643

    set target_resid 540
    set target_resid1 546
    set target_resid2 510
    #The resid you want to analysis
    set r_out 999
    set r_out1 999
    set r_out2 999
    set r_g 999
    #Out_put radius
    set z_out -1
    set z_out1 -1
    set z_out2 -1
    set z_g -1
    set ri -1
    #Out_put z coordinate
    foreach element $hole_list {
	set resid_in_list [lindex $element 3]
        if {$resid_in_list == $target_resid} then {
	    set r_tem [lindex $element 1]
            if { $r_tem < $r_out} then {
		set r_out $r_tem
		set z_out [lindex $element 0]
	    }
	}   elseif {$resid_in_list == $target_resid1} then {
            set r_tem1 [lindex $element 1]
            if { $r_tem1 < $r_out1} then {
                set r_out1 $r_tem1
                set z_out1 [lindex $element 0]
            }
        }   elseif {$resid_in_list == $target_resid2} then {
            set r_tem2 [lindex $element 1]
            if { $r_tem2 < $r_out2} then {
                set r_out2 $r_tem2
                set z_out2 [lindex $element 0]
            }
        }    
        set r_gtem [lindex $element 1]
        if {$r_gtem < $r_g} then {
	    set r_g $r_gtem  
            set z_g [lindex $element 0]
            set ri [lindex $element 3]
        }

    }
    puts $output "$i,$r_out,$z_out"
    puts $output1 "$i,$r_out1,$z_out1"
    puts $output2 "$i,$r_out2,$z_out2"
    puts $outputg "$i,$r_g,$z_g,$ri"
}
close $output
close $output1
close $output2
close $outputg
