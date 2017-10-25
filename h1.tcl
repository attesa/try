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




#Modification for multi-run to reduce the random effect
#Align with s5-s6 first
#5/12/2015

source hole.tcl

set selection [atomselect top protein]
set output [open "1027_p5t303_skip10_671.csv" w]
set output1 [open "1027_p5t303_skip10_679.csv" w]
set output2 [open "1027_p5t303_skip10_643.csv" w]
set output3 [open "1027_p5t303_skip10_lower.csv" w]
set outputg [open "1027_p5t303_skip10_global.csv" w]
set output4 [open "1027_p5t303_skip10_upper.csv" w]

#set output file 
#set output [open "r671" w]
#set output1 [open "r679.dat" w]
#set output2 [open "r643.dat" w]
#set outputg [open "g.dat" w]



set mol_id [$selection molindex]
set num_frames [molinfo $mol_id get numframes]
set target_resid 671
set target_resid1 679
set target_resid2 643
set lower 672



for {set i 0} {$i < $num_frames} {incr i} {
#Loop for all frames
    set z_test -999
    set rt_out -1000
    set rt_out1 -1000
    set rt_out2 -1000
    set rt_out3 -1000
    set rt_out4 -1000
    set rt_g -1000
    set rti -1000
    set rti3 -1000
    set rti4 -1000
    set zt_out -999
    set zt_out1 -999
    set zt_out2 -999
    set zt_out3 -999
    set zt_out4 -999
    set zt_g -999


#    for {set j 0} {$j < 10} {incr j 1} {
    	#loop to reduece random effect
        $selection frame $i
    	puts $i
    	set hole_list [::Hole::runhole $selection]
    	#hole_list contains all the returned value in form of 
    	#[z radius resname resid]
    	#The resid you want to analysis
    	set r_out 999
    	set r_out1 999
    	set r_out2 999
        set r_out3 999
    	set r_g 999
	set r_out4 999
    	#Out_put radius
	set z_out3 -1000
	set z_out4 -1000
    	set z_out -1000
    	set z_out1 -1000
    	set z_out2 -1000
    	set z_g -1000
    	set ri -1000
    	#out_put z coordinate
    	foreach element $hole_list {
	    set resid_in_list [lindex $element 3]
            if {$resid_in_list == $target_resid} then {
	        set r_tem [lindex $element 1]
            	if { $r_tem < $r_out } then {
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
	    if {$resid_in_list > $lower} then {
                set r_tem3 [lindex $element 1]
                if { $r_tem3 < $r_out3} then {
                    set r_out3 $r_tem3
                    set z_out3 [lindex $element 0]
		    set ri3 [lindex $element 3]
                }
    	    } elseif {$resid_in_list < $lower} then {
                set r_tem4 [lindex $element 1]
                if { $r_tem4 < $r_out4} then {
                    set r_out4 $r_tem4
                    set z_out4 [lindex $element 0]
                    set ri4 [lindex $element 3]
                }
            }
 



            set r_gtem [lindex $element 1]
            if {$r_gtem < $r_g} then {
	        set r_g $r_gtem  
                set z_g [lindex $element 0]
                set ri [lindex $element 3]
            }

        }
	#Now compare r of this round with previsous max:
	if {$r_out > $rt_out && $z_out > $z_test }  then {   
	    set rt_out $r_out
	    set zt_out $z_out
	}

        if {$r_out1 > $rt_out1 && $z_out1 > $z_test }  then {
            set rt_out1 $r_out1
            set zt_out1 $z_out1
        }

        if {$r_out2 > $rt_out2 && $z_out2 > $z_test }  then {
            set rt_out2 $r_out2
            set zt_out2 $z_out2
        }

        if {$r_out3 > $rt_out3 && $z_out3 > $z_test }  then {
            set rt_out3 $r_out3
            set zt_out3 $z_out3
	    set rti3 $ri3
        }

        if {$r_out4 > $rt_out4 && $z_out4 > $z_test }  then {
            set rt_out4 $r_out4
            set zt_out4 $z_out4
            set rti4 $ri4
        }



        puts "$r_g,$rt_g"

        if {$r_g > $rt_g && $z_g > $z_test }  then {
            set rt_g $r_g
            set zt_g $z_g
	    set rti $ri
        }



#    }
    #end of j loop






    puts $output4 "$i,$rt_out4,$zt_out4, $rti4"
    puts $output3 "$i,$rt_out3,$zt_out3, $rti3"
    puts $output "$i,$rt_out,$zt_out"
    puts $output1 "$i,$rt_out1,$zt_out1"
    puts $output2 "$i,$rt_out2,$zt_out2"
    puts $outputg "$i,$rt_g,$zt_g,$rti"
}
close $output
close $output1
close $output2
close $outputg
close $output3
close $output4
