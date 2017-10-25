#This script is going to find an optiaml (but not the best) orientation of a protein 
#so that a minmal waterbox can be generated
#
#The principal is to explore all possible orientation with a certain angle step length, by means of Eular angles of z-x-z
#
#Han Wen
#10/20/2015
#




set sel [atomselect top all]

set com [measure center $sel weight mass]
$sel moveby [vecscale -1.0 $com]

set i_op 0
set j_op 0
set k_op 0
set box_size_op 99999999

for {set i 0} {$i < 360} {incr i 1} { 
#rotate around z, third rotation

    set matrix [transaxis z $i]
    $sel move $matrix

	    #measure minmax and calculate box size
	    set minmax_result [measure minmax $sel]
	    set x1 [lindex [lindex $minmax_result 0] 0]
            set y1 [lindex [lindex $minmax_result 0] 1]
            set z1 [lindex [lindex $minmax_result 0] 2]
            set x2 [lindex [lindex $minmax_result 1] 0]
            set y2 [lindex [lindex $minmax_result 1] 1]
            set z2 [lindex [lindex $minmax_result 1] 2]

	    set x [expr {$x2-$x1+30}]
	    set y [expr {$y2-$y1+30}]
            #set z [expr {$z2-$z1+30}]

	    if {$x > $y} {
		set box_size [expr abs($x*$x)]	
	    }

            if {$x < $y} {
                set box_size [expr abs($y*$y)]
            }

	    #set box_size [expr abs($x*$y*$z)]
	    #puts "$i,$j,$k,$box_size"	
            if {$box_size < $box_size_op} {
		set i_op $i 
		set box_size_op $box_size
	    }



    set matrix [transaxis z -$i]
    $sel move $matrix

}


set matrix [transaxis z $i_op]
$sel move $matrix


$sel moveby $com

puts "$i_op, , $box_size_op"

