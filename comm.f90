module comm
 character(3) :: printDB !yes or no
 character(3) :: force2Domains !yes or no
 character(3) :: printMobility !yes or no
 integer :: numVecs
 character(200) :: filename
 character(3) :: doISE,printTCL,printVEC,printCline,printPline,printGline
 integer,parameter :: N_bins=1000
 integer,dimension(1:N_bins) :: SGtable

 ! Compiler such as gfortran and g95 do not understand formats
 ! like "(<N>I4)"
 ! In order to emulate this we need to convert N to string and make format string

 character(20) :: fmt

contains

 function int_to_str(i)
 integer,intent(in) :: i
 character(20) :: int_to_str
  write(int_to_str,"(I10)") i
  int_to_str = trim(int_to_str)
 end function

end module