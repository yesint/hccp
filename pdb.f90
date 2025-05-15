module pdb

type atom_entry
 integer :: atom_num
 character(3) :: atom_type
 character(3) :: res_name
 integer :: res_id
 character(1) :: chain
 integer :: res_num
 real,dimension(1:3) :: xyz
 real :: occupancy,beta
 character(3) :: segname
 real :: mass
end type

!Residue names corresponding to the coding number
character(3),dimension(1:20),parameter :: residue_ids = &
 (/"CYS","MET","PHE","ILE","LEU","VAL","TRP","TYR","ALA","GLY", &
   "THR","SER","GLN","ASN","GLU","ASP","HIS","ARG","LYS","PRO"/)

contains

 subroutine read_pdb(filename,all_atoms,num_atoms,chain)
 character(*),intent(in) :: filename
 character(1) :: chain
 type(atom_entry),dimension(:),pointer :: all_atoms
 integer,intent(out) :: num_atoms
 !===================================
 integer :: ios
 character(6) :: dum1
 character(100) :: dumL
 type(atom_entry) :: temp
 integer :: i,j
 !===================================
  if(chain=="_") chain=" "
  print *,"Reading file ",filename,"..."
  !Scan file from the end to find the number of the last atom
  open(111,file=filename,action="read",iostat=ios,position="append")
!  print *,ios
  if(ios/=0)then
   print *,"File not found!"
   goto 431
  end if

   backspace(111) 
   do while(ios/=-1)
    read(111,"(A6,I5,A)",iostat=ios) dum1,i
    backspace(111) 
    backspace(111)
    !if(trim(dum1)=="ATOM" .or. trim(dum1)=="HETATM")then
    if(trim(dum1)=="ATOM")then
     num_atoms = i
     exit
    end if
   end do
   print *,"Number of atoms:",num_atoms
  !Allocating the array for all atoms
  allocate(all_atoms(1:num_atoms))
  !Reopen the file from the beginning and read the atoms
  rewind(111)
  i = 0
  20 do while(ios/=-1)
    read(111,"(A6,I5,A5,A4,A2,I4,F12.3,F8.3,F8.3,F6.2,F6.2,A10)", iostat=ios,err=20) &
      dum1,temp%atom_num,temp%atom_type,temp%res_name,temp%chain,temp%res_num, &
      temp%xyz(1),temp%xyz(2),temp%xyz(3),temp%occupancy,temp%beta,temp%segname
    !if(trim(dum1)=="ATOM" .or. trim(dum1)=="HETATM")then
    if(trim(dum1)=="ATOM")then
     i = i+1
     !Find residue id
     temp%res_id = 0 !Imply that this is not protein
     do j=1,20
      if(temp%res_name==residue_ids(j))then
       temp%res_id = j
       exit
      end if
     end do
     if(temp%res_id==0) print *,"!!! WARNING: Residue ",temp%res_num,temp%res_name," is not recognized!"
     !Assign the mass of the atoms
     select case(temp%atom_type(1:1))
      case("C")
       temp%mass = 12.0
      case("O")
       temp%mass = 15.0
      case("N")
       temp%mass = 14.0
      case("S")
       temp%mass = 32.0
      case("H")
       temp%mass = 1.0
      case("P")
       temp%mass = 31.0
     end select
     !Restrict to chain
     if(temp%chain==chain .or. chain==" ")then
      all_atoms(i) = temp
     else
      i = i-1
     end if
     !print *,i,all_atoms(i)%atom_num," ",all_atoms(i)%atom_type,all_atoms(i)%res_num,all_atoms(i)%xyz
    end if
    if(i==num_atoms) exit
  end do
  close(111)
  !Reading done
  print *,"Finished with file ",filename
431 end subroutine

 subroutine write_pdb(filename,atoms)
 character(*),intent(in) :: filename
 type(atom_entry),dimension(:),pointer :: atoms
 !---------------------------
 integer :: N,i
 !---------------------------
  N = size(atoms)
  open(111,file=trim(filename),action="write")
  !Write header
  !
  !Write atoms
  do i=1,N
   write(111,"(A6,I5,A5,A4,A2,I4,F12.3,F8.3,F8.3,F6.2,F6.2,A10)") &
      "ATOM  ",atoms(i)%atom_num,adjustr(atoms(i)%atom_type),adjustr(atoms(i)%res_name),adjustr(atoms(i)%chain),atoms(i)%res_num, &
      atoms(i)%xyz(1),atoms(i)%xyz(2),atoms(i)%xyz(3),atoms(i)%occupancy,atoms(i)%beta,adjustr(atoms(i)%segname)
  end do
  close(111)
 end subroutine

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 
 subroutine select_atoms(all_atoms,selected_atoms,selection,num_selected)
 type(atom_entry),dimension(:),pointer :: all_atoms,selected_atoms
 character(*),intent(in) :: selection
 integer,intent(out) :: num_selected
 !===================================
 integer :: i,j,num_atoms,cur_res,ca_found
 real,dimension(1:3) :: cm
 real :: res_mass,res_beta,res_sz
 character(10) :: ch
 !===================================
 num_atoms = size(all_atoms)
 print *,"Selecting ",selection," among ",num_atoms," atoms..."
  !Extracting information about selected atoms
  num_selected = 0
  if(trim(selection)=="SCM" .or. trim(selection)=="SCG")then
   cur_res = 0
   i=1
   num_selected = 0
   mm: do
    ca_found = 0
    do while(all_atoms(i)%res_num==cur_res)
     !print *,"  --",i,all_atoms(i)%res_name,all_atoms(i)%res_num,all_atoms(i)%atom_type
     if( all_atoms(i)%atom_type(1:2)=="CA" )then
      ca_found = 1
     end if
     i = i+1
     if(i>num_atoms) goto 565
    end do
    cur_res = all_atoms(i)%res_num
565 if(i/=1)then
     if(all_atoms(i-1)%res_id/=0 .and. ca_found==1) then
      !print *,">>>",all_atoms(i-1)%res_num,all_atoms(i-1)%res_name,i-1
      num_selected = num_selected + 1
     end if
     if(i>num_atoms) exit mm
    end if
    
   end do mm

  else
   !print *,"Ordinary selection by atom type ",selection
   do i=1,num_atoms
    if(trim(all_atoms(i)%atom_type)==trim(selection) .and. all_atoms(i)%res_id/=0 )then 
     num_selected = num_selected+1
     !print *,num_selected,all_atoms(i)%atom_num,all_atoms(i)%atom_type,all_atoms(i)%res_num
    end if
   end do
  end if
  print *,"Found ",num_selected," ",selection," atoms"

  !Allocating array for selected atoms
  !print *,"!!!",num_selected
  allocate(selected_atoms(1:num_selected))
  !print *,"!!!",num_selected
  !Fill the array
  if(trim(selection)=="SCM" .or. trim(selection)=="SCG")then
   !Special case: SCM - Side chain center of masses
   cur_res = 0
   i = 1
   j = 1
   do !Cycle over all residues
    
    !print *,cur_res
    !print *,"!!!",num_selected
    !
    if(j>num_selected)then
     !print *,"Exting via j ",j,num_selected
     exit
    end if
    cm = 0
    res_mass = 0
    res_beta = 0
    res_sz = 0
    ca_found = 0
    do while(all_atoms(i)%res_num==cur_res .and. i<=num_atoms)

     ch = all_atoms(i)%atom_type
     !!print *,"Atom ",i," ",all_atoms(i)%res_name,"(",all_atoms(i)%res_num,") ",ch

     if(all_atoms(i)%res_name=="GLY")then
      !!print *,"  GLY residue detected"
      i = i+1
      if( ch(1:2)=="CA" )then
        !!print *,"  Setting coordinates for CA since no sidechain is present"
        ca_found = 1
        res_sz = res_sz + 1
        if(trim(selection)=="SCG")then
         cm(1:3) = cm(1:3) + all_atoms(i)%xyz(1:3)
        else
         cm(1:3) = cm(1:3) + all_atoms(i)%xyz(1:3) * all_atoms(i)%mass
        end if
        res_mass = res_mass + all_atoms(i)%mass
        res_beta = res_beta + all_atoms(i)%beta
      end if
      goto 545
     end if
     
     !Only non-hydrogens and not backbone
     if( ch(1:2)=="CA" ) ca_found = 1
     if( ch(1:1)/="H" .and. ch(1:2)/="O " .and. ch(1:2)/="C " .and. ch(1:2)/="N " .and. ch(1:2)/="CA")then 
      !!print *,"  Adding coordinates for this atom"
      res_sz = res_sz + 1
      !Sum up coordinates
      if(trim(selection)=="SCG")then
       cm(1:3) = cm(1:3) + all_atoms(i)%xyz(1:3)
      else
       cm(1:3) = cm(1:3) + all_atoms(i)%xyz(1:3) * all_atoms(i)%mass
      end if
      !Sum up masses
      res_mass = res_mass + all_atoms(i)%mass
      res_beta = res_beta + all_atoms(i)%beta
      i = i+1
     else
      i = i+1
     end if
545 end do
    
    !print *,"End of residue"
    
    !End of the residue
    if(i/=1 .and. ca_found==1)then
     selected_atoms(j) = all_atoms(i-1)
     selected_atoms(j)%mass = res_mass
     selected_atoms(j)%beta = res_beta/res_sz
     if(trim(selection)=="SCG")then
      selected_atoms(j)%xyz = cm/res_sz
     else
      selected_atoms(j)%xyz = cm/res_mass
     end if
     selected_atoms(j)%atom_type = "CA"
     !!print *,"coor:",selected_atoms(j)%xyz
     !!print *," Residue ",selected_atoms(j)%res_name,selected_atoms(j)%res_num,j
     j = j+1
     !!print *,j
    end if
    !Set next res_id
    if(i<=num_atoms)then
     cur_res = all_atoms(i)%res_num
    else
     !That's all! Exit.
     exit
    end if
   end do
  else
   !Ordinary selection
   j=0
   do i=1,num_atoms
    if(trim(all_atoms(i)%atom_type)==trim(selection)  .and. all_atoms(i)%res_id/=0)then
     j=j+1
     selected_atoms(j) = all_atoms(i)
     !print *,j,selected_atoms(j)%res_name,selected_atoms(j)%res_num
    end if
   end do
  end if
  !print *,"Seletion done."
 end subroutine

 subroutine apply_from_ref2(atoms,ref,vectors,mode1,mode2,s1,s2)
 type(atom_entry),dimension(:),pointer :: atoms,ref
 real,dimension(:,:),pointer :: vectors
 integer :: vec_num
 real :: s1,s2
 integer :: mode1,mode2
 !-------------------------------------------------
 integer :: i,N
 !-------------------------------------------------
  N = size(atoms)
  do i=1,N
   atoms(i)%xyz(1) = ref(i)%xyz(1) + s1*vectors(3*(i-1)+1,mode1) + s2*vectors(3*(i-1)+1,mode2)
   atoms(i)%xyz(2) = ref(i)%xyz(2) + s1*vectors(3*(i-1)+2,mode1) + s2*vectors(3*(i-1)+2,mode2)
   atoms(i)%xyz(3) = ref(i)%xyz(3) + s1*vectors(3*(i-1)+3,mode1) + s2*vectors(3*(i-1)+3,mode2)
  end do
 end subroutine

 subroutine apply_from_ref(atoms,ref,vectors,mode1,s1)
 type(atom_entry),dimension(:),pointer :: atoms,ref
 real,dimension(:,:),pointer :: vectors
 integer :: vec_num
 real :: s1
 integer :: mode1
 !-------------------------------------------------
 integer :: i,N
 !-------------------------------------------------
  N = size(atoms)
  do i=1,N
   atoms(i)%xyz(1) = ref(i)%xyz(1) + s1*vectors(3*(i-1)+1,mode1)
   atoms(i)%xyz(2) = ref(i)%xyz(2) + s1*vectors(3*(i-1)+2,mode1)
   atoms(i)%xyz(3) = ref(i)%xyz(3) + s1*vectors(3*(i-1)+3,mode1)
  end do
 end subroutine

 subroutine apply_vec(atoms,ref,vectors,vec)
 type(atom_entry),dimension(:),pointer :: atoms,ref
 real,dimension(:,:),pointer :: vectors
 real,dimension(:),pointer :: vec
 !-------------------------------------------------
 integer :: i,N,N_vec
 real :: x,y,z
 !-------------------------------------------------
  N_vec = size(vec)
  N = size(ref)
  do i=1,N !Atoms
   x=0;  y=0;  z=0
   !Accumulate shifts from all vectors
   do j=1,N_vec  !Modes
    x = x+vec(j)*vectors(3*(i-1)+1,j)
    y = y+vec(j)*vectors(3*(i-1)+2,j)
    z = z+vec(j)*vectors(3*(i-1)+3,j)
   end do
   !Apply shifts
   atoms(i)%xyz(1) = ref(i)%xyz(1) + x
   atoms(i)%xyz(2) = ref(i)%xyz(2) + y
   atoms(i)%xyz(3) = ref(i)%xyz(3) + z
  end do
 end subroutine


 subroutine apply_eigenvector(atoms,vectors,vec_num,scale)
 type(atom_entry),dimension(:),pointer :: atoms
 real,dimension(:,:),pointer :: vectors
 integer :: vec_num
 real :: scale
 !-------------------------------------------------
 integer :: i,N
 !-------------------------------------------------
  N = size(atoms)
  do i=1,N
   atoms(i)%xyz(1) = atoms(i)%xyz(1) + scale*vectors(3*(i-1)+1,3*N-5-vec_num)
   atoms(i)%xyz(2) = atoms(i)%xyz(2) + scale*vectors(3*(i-1)+2,3*N-5-vec_num)
   atoms(i)%xyz(3) = atoms(i)%xyz(3) + scale*vectors(3*(i-1)+3,3*N-5-vec_num)
  end do
 end subroutine

 subroutine apply_by_residue(all_atoms,selected_atoms,vectors,vec_num,scale)
 type(atom_entry),dimension(:),pointer :: all_atoms,selected_atoms
 real,dimension(:,:),pointer :: vectors
 integer :: vec_num
 real :: scale
 !-------------------------------------------------
 integer :: i,j,N,Nall
 !-------------------------------------------------
  N = size(selected_atoms) !N - number of residues
  Nall = size(all_atoms) !number of atoms
  do i=1,N !Cycle over residues
   do j=1,Nall 
    !Find atoms from residue i and move them
    if(all_atoms(j)%res_num==selected_atoms(i)%res_num)then
     all_atoms(j)%xyz(1) = all_atoms(j)%xyz(1) + scale*vectors(3*(i-1)+1,3*N-5-vec_num)
     all_atoms(j)%xyz(2) = all_atoms(j)%xyz(2) + scale*vectors(3*(i-1)+2,3*N-5-vec_num)
     all_atoms(j)%xyz(3) = all_atoms(j)%xyz(3) + scale*vectors(3*(i-1)+3,3*N-5-vec_num)
    end if
   end do
  end do
 end subroutine

end module