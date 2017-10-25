#Simple script to iteratively measure com 
#Han Wen
#9/29/2016
#

set selection [atomselect top "resid 405 to 535 and name CA and chain A"]
set mol_id [$selection molindex]
set com_out [open "com_test.csv" w]
set num_frames [molinfo $mol_id get numframes]
puts $num_frames
for {set i 0} {$i < $num_frames} {incr i} {
    set com [measure center $selection weight mass]
    set x [lindex $com 0]	
    set y [lindex $com 0]
    set z [lindex $com 0]
    puts $com_out "$i,$x,$y,$z"
}
close $com_out
