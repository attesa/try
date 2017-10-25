#measure 3j5p twist angle, just recording down x,y coordinate and the angle
#hanwen
#



set A [atomselect top "resid 111 to 280 and index 0 to 9533"]
set B [atomselect top "resid 111 to 280 and index 9534 to 19067"]
set C [atomselect top "resid 111 to 280 and index 19068 to 28601"]
set D [atomselect top "resid 111 to 280 and index 28602 to 38135"]


set angle_outA [open "sp1_twist_A.csv" w]
set angle_outB [open "sp1_twist_B.csv" w]
set angle_outC [open "sp1_twist_C.csv" w]
set angle_outD [open "sp1_twist_D.csv" w]


set mol_id [$A molindex]
set num_frames [molinfo $mol_id get numframes]
puts $num_frames
for {set i 0} {$i < $num_frames} {incr i} {
    $A frame $i
    $B frame $i
    $C frame $i
    $D frame $i
    set CA [measure center $A]
    set CB [measure center $B]
    set CC [measure center $C]
    set CD [measure center $D]
    set CA_x [lindex $CA 0]
    set CA_y [lindex $CA 1]
    set CB_x [lindex $CB 0]
    set CB_y [lindex $CB 1]
    set CC_x [lindex $CC 0]
    set CC_y [lindex $CC 1]
    set CD_x [lindex $CD 0]
    set CD_y [lindex $CD 1]
    set tanA [expr {$CA_y/$CA_x}]
    set tanB [expr {$CB_y/$CB_x}]
    set tanC [expr {$CC_y/$CC_x}]
    set tanD [expr {$CD_y/$CD_x}]
    set an_A [expr {180 * atan($tanA) / 3.1415926}]
    set an_B [expr {180 * atan($tanB) / 3.1415926}]
    set an_C [expr {180 * atan($tanC) / 3.1415926}]
    set an_D [expr {180 * atan($tanD) / 3.1415926}]

    puts $angle_outA "$i, $an_A"
    puts $angle_outB "$i, $an_B"
    puts $angle_outC "$i, $an_C"
    puts $angle_outD "$i, $an_D"
}
close $angle_outA
close $angle_outB
close $angle_outC
close $angle_outD




