#Simple script to iteratively measure com 
#Han Wen
#9/29/2016
#
#further modified to measure the com distance between two domain
#
#
#
#10/2/2016


#set selection1 [atomselect top "resid 393 to 534 757 to 787 and name CA and chain C"]
#set selection2 [atomselect top "resid 535 to 544 662 to 756 788 to 798 and name CA and chain C"]


#set selection1 [atomselect top "resid 400 to 529 760 to 790 and name CA and chain B"]
#set selection2 [atomselect top "resid 530 to 541 660 to 759 791 to 800 and name CA and chain B"]


#set selection1 [atomselect top "resid 484 485 and chain A"]
#set selection2 [atomselect top "resid 688 689 and chain A"]

set selection1 [atomselect top "resid 405 406 407 and chain A"]
set selection2 [atomselect top "resid 714 715 and chain A"]

#set selection1 [atomselect top "resid 485 486 and chain B"]
#set selection2 [atomselect top "resid 689 690 and chain B"]

#set selection1 [atomselect top "resid 413 414 415  and chain B"]
#set selection2 [atomselect top "resid 713 714 and chain B"]



set mol_id [$selection1 molindex]
set com_out [open "com_dist.csv" w]
set num_frames [molinfo $mol_id get numframes]
puts $num_frames
for {set i 0} {$i < $num_frames} {incr i} {
    $selection1 frame $i
    $selection2 frame $i

    set com1 [measure center $selection1 weight mass]
    set x1 [lindex $com1 0]	
    set y1 [lindex $com1 1]
    set z1 [lindex $com1 2]
    set com2 [measure center $selection2 weight mass]
    set x2 [lindex $com2 0]
    set y2 [lindex $com2 1]
    set z2 [lindex $com2 2]


    set dist [expr {sqrt(($x1-$x2) * ($x1-$x2)  + ($y1-$y2) * ($y1-$y2)  +($z1-$z2) * ($z1-$z2))}]

    puts $com_out "$i $dist"
}
close $com_out
