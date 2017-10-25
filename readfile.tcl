set file [open "try.tcl" r]
set number 0
while { [gets $file line] >=0} {
    incr number
    puts "[gets $file line]"	
}
close $file
puts "Number of lines: $number" 
