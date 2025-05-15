module input_parser

use input_parser_util

implicit none
 
 interface define_keyword
 
  subroutine define_keyword_integer(name,group,data_type,variable,default_value)
  use input_parser_util
  implicit none
  character(*) :: name,group,data_type
  integer,intent(out) :: variable
  integer,optional :: default_value
  end subroutine
 
  subroutine define_keyword_real(name,group,data_type,variable,default_value)
  use input_parser_util
  implicit none
  character(*) :: name,group,data_type
  real,intent(out) :: variable
  real,optional :: default_value
  end subroutine
 
  subroutine define_keyword_string(name,group,data_type,variable,default_value,options)
  use input_parser_util
  implicit none
  character(*) :: name,group,data_type
  character(*),intent(out) :: variable
  character(*),optional :: default_value
  character(*),dimension(:),optional,intent(in) :: options
  end subroutine

 end interface

contains

 subroutine parse_file(file_name)
 character(*) :: file_name
 !-------------------------------
 integer :: i,j,L,k,ok,ios,keyword_index,line
 character(80) :: str,word,cur_group
 !-------------------------------
  print *,"(I) INPUT PARSER STARTED "

  cur_group=""
  errors = 0
  num_of_keywords = 0
  line = 0
  !See if don't we have stdin as parameter
  if(trim(file_name)/="<STDIN>")then
   !Open file from disk
   open(111,file=trim(file_name),action="read",iostat=ios)
  end if

  do 
   !Read the line
   if(trim(file_name)=="<STDIN>")then
    read(*,"(A)",iostat=ios) str
   else
    read(111,"(A)",iostat=ios) str
   end if

   if(ios == -1) exit !End of file

   line = line + 1

   !Parse the string
   str = trim(adjustl(str))
   L = len(str) !Length of string
   !If first character is # or the line is empty skip this line
   if(str(1:1)/="#" .and. trim(str)/="")then
    !Find the first word in the string
    k = index(str," ")
    word = str(1:k-1) 
    !print *,"first word> ",trim(word)
    !If this word is the name of the group?
    
    if( trim(str(k+1:L))=="" .and. trim(word)/="end" )then
     !Yes, this is group name
     if(trim(cur_group)=="")then
      cur_group = word
      !print *,"> Group ",trim(cur_group)," is started."
     end if
    else
     !No, this is not a group name. 
     !Is it end directive?
     if(trim(word)=="end")then
      !print *,"> Group ",trim(cur_group)," is finished."
      cur_group = ""
     else
      !This can only be a keyword 
      !Processing this keyword
      !print *,"> Processing keyword ",trim(word)
      j = get_index_by_name(trim(word))
      if(j/=0)then
       !This keyword was already present above
       print "(A,A,A,I4,A)"," (E) Keyword ",trim(word)," was already defined above :(  [",line,"]"
       errors = errors + 1
      else
       num_of_keywords = num_of_keywords + 1
       keywords(num_of_keywords)%name = trim(word)
       keywords(num_of_keywords)%group_name = trim(cur_group)
       keywords(num_of_keywords)%val = trim(adjustl(str(k+1:L)))
      end if
     end if
    end if !End of group test

   end if !End of comment test

  end do !End of reading file

  if(trim(file_name)/="<STDIN>") close(111)

 end subroutine


 subroutine validate_keywords()
 logical :: exit_on_error
 integer :: i
  !Cycle over all keywords and see if all of them are defined
  do i=1,num_of_keywords
   if(keywords(i)%ok .eqv. .false.)then
    print *,"(W) Keyword ",trim(keywords(i)%name)," is not recognized and ignored"
   end if
  end do

  if(errors>0)then
   print "(A,I4,A)", " (E) There are ",errors," errors, further processing can be unpredictable :("
  else
   print *, "(I) Input parsed with no errors :)"
  end if
  print *,"(I) INPUT PARSER FINISHED "

 end subroutine

end module

 !Assigning value routines must be outside of the module's scope

  subroutine define_keyword_integer(name,group,data_type,variable,default_value)
  use input_parser_util
  implicit none
  character(*) :: name,group,data_type
  integer,intent(out) :: variable  
  integer,optional :: default_value
  integer :: ios
  logical :: found
   !Run default action
   call define_keyword_general(name,group,found)
   if(found .eqv. .true.)then
    read(keywords(get_index_by_name(name))%val,*,iostat=ios) variable
    if(ios/=0)then
     !Failure :(
     variable = 0
     print *,"(E) Value '",trim(keywords(get_index_by_name(name))%val),"' is not valid integer for keyword '",name,"' :("
     errors = errors + 1
    end if
   else
    if(.not. present(default_value))then
     print *,"(E) Keyword ",trim(name)," is not found!"
     errors = errors + 1
     variable = 0
    else
     print *,"(I) Keyword ",trim(name)," is not found. Using default ",default_value
     variable = default_value
    end if
   end if
  end subroutine

  subroutine define_keyword_real(name,group,data_type,variable,default_value)
  use input_parser_util
  implicit none
  character(*) :: name,group,data_type
  real,optional :: default_value
  real,intent(out) :: variable  
  integer :: ios
  logical :: found
   !Run default action
   call define_keyword_general(name,group,found)
   if(found .eqv. .true.)then
    read(keywords(get_index_by_name(name))%val,*,iostat=ios) variable
    if(ios/=0)then
     !Failure :(
     variable = 0
     print *,"(E) Value '",trim(keywords(get_index_by_name(name))%val),"' is not valid float for keyword '",name,"' :("
     errors = errors + 1
    end if
   else
    if(.not. present(default_value))then
     print *,"(E) Keyword ",trim(name)," is not found!"
     errors = errors + 1
     variable = 0.0
    else
     print *,"(I) Keyword ",trim(name)," is not found. Using default ",default_value
     variable = default_value
    end if
   end if
  end subroutine

  subroutine define_keyword_string(name,group,data_type,variable,default_value,options)
  use input_parser_util
  implicit none
  character(*) :: name,group,data_type
  character(*),intent(out) :: variable  
  character(*),optional :: default_value
  character(*),dimension(:),optional,intent(in) :: options
  integer :: ios,i,j
  logical :: found
   !Run default action
   call define_keyword_general(name,group,found)
   if(found .eqv. .true.)then
    variable = trim(keywords(get_index_by_name(name))%val)
    !If options are given see if the value is valid
    if(present(options))then
     j = 0
     do i=1,size(options)
      if(trim(variable)==trim(options(i)))then 
       j=i
       exit
      end if
     end do
     !If within options do...
     if(j==0)then
      print *,"(E) The value '",trim(variable),"' is not valid for keyword '",trim(keywords(get_index_by_name(name))%name),"'"
      print *,"    Allowed values are:"
      do i=1,size(options)
       print *,"    '",trim(options(i)),"'"
      end do
     end if
    end if
   else
    if(.not. present(default_value))then
     print *,"(E) Keyword ",trim(name)," is not found!"
     errors = errors + 1
     variable = ""
    else
     print *,"(I) Keyword ",trim(name)," is not found. Using default ",default_value
     variable = default_value
    end if
   end if
  end subroutine

