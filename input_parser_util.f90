!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  THIS MODULE IS A PART OF THE INPUT_PARSER - THE FACILITY, WHICH
!  ALLOWS TO PARSE SIMPLE INPUT FILES.
!
!  THIS FILE SHOULD NEVER BE USEd DIRECTLY IN MAIN PROGRAM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  CONTAINS UTILITES FOR THE INPUT PARSER AND DATA STRUCTURES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  (C) 2005, Semen Yesylevskyy
!  V. 2.0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module input_parser_util

 type keyword_description
  character(80) :: group_name
  character(80) :: name
  character(80) :: data_type
  character(200) :: val
  logical :: ok = .false.
 end type

 type(keyword_description),dimension(1:100),save :: keywords
 integer :: num_of_keywords=0
 integer :: errors

contains

  function get_index_by_name(keyword)
  character(*) :: keyword
  integer :: get_index_by_name
  integer :: i
   get_index_by_name = 0
   do i=1,num_of_keywords
    if(trim(keyword)==trim(keywords(i)%name))then
     get_index_by_name = i
     exit
    end if
   end do
  end function

  !Checks the name of the keyword, group and other things
  subroutine define_keyword_general(name,group,found)
  character(*) :: name,group
  logical,intent(out) :: found
  integer :: i,j
   !Find this keyword
   !print *,"> defining ",name," ",group
   found = .true.
   i = get_index_by_name(trim(name))
   if(i==0)then
    !Keyword not found
    found = .false.
   else
    !Keyword recognized
    keywords(i)%ok = .true.
    !Check if group is valid
    !If group = "" then this keyword can be everywhere
    if(trim(keywords(i)%group_name)/=trim(group) .and. trim(group)/="" )then
     !Not valid
     print *,"(E) Keyword ",trim(name)," is only allowed in group ",trim(keywords(i)%group_name),", not in ",trim(group)," :("
     errors = errors + 1
    end if
   end if
  end subroutine

end module