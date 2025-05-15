module functions
use pdb
use eigen
use comm
implicit none

!Conversion factor to convert DFIRE values to kT values
!real,parameter :: dfire_factor=0.0157 * 4.187e3 / ( 6.02e23 * 300.0 * 1.380658e-23 )
real,parameter :: dfire_factor=1.0
!                              ^        ^           ^         ^       ^
!                              eta      converts   Avagadro   ---------
!                      (to Kcal/Mole)   to J/Mole   (to J)       kT
contains

subroutine do_gnm(N,rc,atoms,values,vectors)
integer :: N
real :: rc
type(atom_entry),dimension(:),pointer :: atoms
real,dimension(:),pointer :: values
real,dimension(:,:),pointer :: vectors
!=============================================
integer :: i,j,k
real :: prom,time1,time2
real,dimension(:),allocatable :: vec2
!real,dimension(:,:),allocatable :: g,g1
!=============================================
allocate(vec2(1:N))
!allocate(g(1:N,1:N),g1(1:N,1:N))
print *,"Constructing GNM Kirkgoff matrix. Cut-off:",rc," Size:",N
vectors = 0
do i=1,N-1
 do j=i+1,N
  prom = sum( ( atoms(i)%xyz(1:3)-atoms(j)%xyz(1:3) )**2 )
  if(sqrt(prom)<=rc)then
   vectors(i,j)= -1
   vectors(j,i)= -1
  end if
 end do
end do

 if(printGline=="yes")then
  fmt = "(A,"//int_to_str(N*N)//"I1)"
  open(111,file="g_string.dat",action="write",access="append")
   write(111,fmt) trim(filename)//" ",int(abs(vectors))
  close(111)
 end if

!Compute diagonal elements
do i=1,N  
 vectors(i,i)=-sum(vectors(i,:))
end do
call cpu_time(time1)

 !print *,"Writing g-matrix..."
 !open(111,file="p.dat",action="write")
 ! do i=1,N
 !  write(111,"(<N>F10.5)") vectors(i,1:N)
 ! end do
 !close(111)

 print *,"Transforming matrix..."
 call tred2(vectors,values,vec2)
 print *,"Computing eigenvectors..."
 call tqli(values,vec2,vectors)


! g = vectors
!    open(111,file="e:\md\gnm\code\g1.dat",action="read")
!     read(111,"(<N>F12.5)") g1
!    close(111)
!    open(111,file="g.dat",action="write")
!     write(111,"(<N>F12.5)") g - g1
!    close(111)

 !print *,"Computing eigenvectors..."
 !call jacobi(g,values,vectors,k)
call cpu_time(time2)
print "(A,F10.5)","Time of finding eigenvectors: ",time2-time1

print *,"Sorting eigenvectors..."
call eigsrt(values,vectors)

!deallocate(vec2,g)
deallocate(vec2)
end subroutine

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine do_anm(N,rc,atoms,values,vectors,weight)
integer :: N
real :: rc
type(atom_entry),dimension(:),pointer :: atoms
real,dimension(:),pointer :: values
real,dimension(:,:),pointer :: vectors
logical,optional,intent(in) :: weight
!=============================================
integer :: i,j,k,i1,j1
real :: prom,time1,time2
real,dimension(:),allocatable :: vec2
real :: m
!=============================================
allocate(vec2(1:3*N))
print *,"Calculating ANM Hessian matrix. Cut-off:",rc," Size:",3*N
if(present(weight)) print *,"Hessian is mass-weighed."
vectors = 0
do i=1,N-1 !Counts super-elements
 do j=i+1,N
  !print *,"Super-element",i,j
  prom = sum( ( atoms(i)%xyz(1:3)-atoms(j)%xyz(1:3) )**2 )
  if(sqrt(prom)<=rc)then
   !If weighting is requested calculate weigthing term
   if(present(weight))then
    m = 1.0/sqrt(atoms(i)%mass * atoms(j)%mass)
   else
    m = 1.0
   end if
   !Fill individual elements of the superelement (i,j)
   do i1=1,3
    do j1=1,3
     vectors(3*(i-1)+i1,3*(j-1)+j1) = -m*( atoms(j)%xyz(i1)-atoms(i)%xyz(i1) )*( atoms(j)%xyz(j1)-atoms(i)%xyz(j1) )/prom;
     vectors(3*(j-1)+j1,3*(i-1)+i1) = vectors(3*(i-1)+i1,3*(j-1)+j1)
    end do
   end do
  end if
 end do
end do
!Compute diagonal elements
do i=1,N !Counts superelements
 !Fill elements of the superelement
 do i1=1,3
  do j1=1,3
   prom = 0
   do k=1,N  
    if(k/=i) prom= prom - vectors(3*(i-1)+i1,3*(k-1)+j1)
   end do
   vectors(3*(i-1)+i1,3*(i-1)+j1)= prom
  end do
 end do
end do
!Starting diagonalization
call cpu_time(time1)
 print *,"Transforming matrix..."
 call tred2(vectors,values,vec2)
 print *,"Computing eigenvectors..."
 call tqli(values,vec2,vectors)
call cpu_time(time2)
print "(A,F10.5)","Time of finding eigenvectors: ",time2-time1

print *,"Sorting eigenvectors..."
call eigsrt(values,vectors)

deallocate(vec2)
end subroutine

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine read_dfire(filename,pot)
character(*),intent(in) :: filename
real,dimension(:,:,:),pointer :: pot
!=================================
integer,parameter :: bin_num=20
integer :: ios,i,j,k
character(3) :: ch1,ch2
real :: val
!=================================
 !Open the file and read hash and potential
 open(111,file=filename,action="read",iostat=ios)
  do while(ios/=-1)
   read(111,*,iostat=ios) ch1,ch2,val,i,j,k
   pot(j+1,k+1,i) = val
   pot(k+1,j+1,i) = val
  end do
 close(111)
end subroutine


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function dfire_energy(atoms,pot)
type(atom_entry),dimension(:),pointer :: atoms
real,dimension(:,:,:),pointer :: pot
real :: dfire_energy
!===================================
integer :: N,i,j,k
real :: r
!===================================
 N = size(atoms)
 dfire_energy = 0
 do i=1,N-1
  do j=i+1,N
   !Calculate the distance between the atoms
   r = sqrt( sum((atoms(i)%xyz-atoms(j)%xyz)**2) )
   !Find needed bin
   if(r<2.0)then
    k=1
   else if(r>=2.0 .and. r<8.0)then
    k = int((r-2.0)/0.5)+2
   else if(r>=8.0 .and. r<15.0)then
    k = int(r-8.0)+14
   else
    k = 20
   end if

   !print *,atoms(i)%res_id,atoms(i)%res_name,atoms(j)%res_id,atoms(j)%res_name
   !print *,i,j,k,"--",r,N
   !print *,atoms(i)%res_id,atoms(j)%res_id
   !!
     if(atoms(i)%res_id==0) print *,"Zero resid at",atoms(i)%res_num,atoms(i)%res_name
     if(atoms(j)%res_id==0) print *,"Zero resid at",atoms(j)%res_num,atoms(j)%res_name
   !!
   dfire_energy = dfire_energy + pot(atoms(i)%res_id,atoms(j)%res_id,k)
  end do
 end do
 dfire_energy = dfire_energy * dfire_factor !Returns energy in kT
end function


function dfire_energy_index(atoms,cl,pot)
integer,dimension(:),pointer :: cl
type(atom_entry),dimension(:),pointer :: atoms
real,dimension(:,:,:),pointer :: pot
real :: dfire_energy_index
!===================================
integer :: N,i,j,k
real :: r
!===================================
 N = count(cl>0)
 dfire_energy_index = 0
 do i=1,N-1
  do j=i+1,N
   !Calculate the distance between the atoms
   r = sqrt( sum((atoms(cl(i))%xyz-atoms(cl(j))%xyz)**2) )
   !Find needed bin
   if(r<2.0)then
    k=1
   else if(r>=2.0 .and. r<8.0)then
    k = int((r-2.0)/0.5)+2
   else if(r>=8.0 .and. r<15.0)then
    k = int(r-8.0)+14
   else
    k = 20
   end if

   !print *,atoms(i)%res_id,atoms(i)%res_name,atoms(j)%res_id,atoms(j)%res_name
   !print *,i,j,k,"--",r,N
   !print *,atoms(i)%res_id,atoms(j)%res_id
   !!
   !  if(atoms(i)%res_id==0) print *,"Zero resid at",atoms(i)%res_num,atoms(i)%res_name
   !  if(atoms(j)%res_id==0) print *,"Zero resid at",atoms(j)%res_num,atoms(j)%res_name
   !!
   dfire_energy_index = dfire_energy_index + pot(atoms(cl(i))%res_id,atoms(cl(j))%res_id,k)
  end do
 end do
 dfire_energy_index = dfire_energy_index * dfire_factor !Returns energy in kT
end function

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function corr(v1,v2)
real,dimension(:) :: v1,v2
real :: corr
!----------------------------------
integer :: N,i,j,m1,m2,s1,s2
!----------------------------------
 N = size(v1)
 m1 = sum(v1)/(N-1)
 m2 = sum(v2)/(N-1)
 s1 = sqrt(sum( (v1-m1)**2 ))
 s2 = sqrt(sum( (v2-m2)**2 ))
 corr = sum( (v1-m1)*(v2-m2)  )/(s1*s2)
end function

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function rmsd(conf1,conf2)
real,dimension(:,:) :: conf1,conf2
real :: rmsd
integer :: N,i
!-----------------------------------
 N = size(conf1(:,1))
 rmsd = 0
 do i=1,N
  rmsd = rmsd + sum((conf1(i,1:3)-conf2(i,1:3))**2)
 end do
 rmsd = sqrt(rmsd/N)
end function
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine make_corr_matrix_ANM(vectors,values,c0,normalize)
real,dimension(:),pointer :: values
real,dimension(:,:),pointer :: vectors,c0
logical,optional,intent(in) :: normalize 
!----------------------------------------
integer :: i,j,k,N
real,dimension(:,:),allocatable :: temp
!----------------------------------------
 N = int(size(c0(1,:))) !c0 is of size N, vectors - 3*N
 allocate(temp(1:3*N,1:3*N))

 print *,"Calculating correlations..."
 do k=1,3*N-6 !For all eigenvectors
! do k=3*N-6,3*N-6
  do i=1,3*N
   do j=i,3*N
    temp(i,j) = temp(i,j) + vectors(i,k)*vectors(j,k)/values(k) 
   end do
  end do
 end do
  
  do i=1,3*N-1
   do j=i,3*N
    temp(j,i) = temp(i,j)
   end do
  end do

  print *,"Extracting diagonal elements..."
  do i=1,N !Cycle over superelements
   do j=1,N
    c0(i,j) = temp(3*(i-1)+1,3*(j-1)+1) + temp(3*(i-1)+2,3*(j-1)+2) + temp(3*(i-1)+3,3*(j-1)+3)

     !c0(i,j) =temp(3*(i-1)+1,3*(j-1)+1)*temp(3*(i-1)+2,3*(j-1)+2)*temp(3*(i-1)+3,3*(j-1)+3)/ &
     !         abs(temp(3*(i-1)+1,3*(j-1)+1)*temp(3*(i-1)+2,3*(j-1)+2)*temp(3*(i-1)+3,3*(j-1)+3))* &  
     !         (abs(temp(3*(i-1)+1,3*(j-1)+1)) + &
     !          abs(temp(3*(i-1)+2,3*(j-1)+2)) + &
     !          abs(temp(3*(i-1)+3,3*(j-1)+3)))
     !c0(i,j) = temp(3*(i-1)+1,3*(j-1)+1)
   end do
  end do

 if(present(normalize))then
  if(normalize .eqv. .true.)then
   do i=1,N
    do j=1,N
     c0(i,j) = c0(i,j)/sqrt(c0(i,i)*c0(j,j)) 
    end do
   end do
  end if
 end if

deallocate(temp)
end subroutine
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine make_corr_matrix_GNM(vectors,values,c0,normalize)
real,dimension(:),pointer :: values
real,dimension(:,:),pointer :: vectors,c0
logical,optional,intent(in) :: normalize 
!----------------------------------------
integer :: i,j,k,N
!----------------------------------------
 N = int(size(values))
 print *,"Calculating correlations..."
 do k=1,N-1 !For all eigenvectors
  do i=1,N
   do j=i,N
    c0(i,j) = c0(i,j) + vectors(i,k)*vectors(j,k)/values(k) 
   end do
  end do
 end do
  
  do i=1,N-1
   do j=i,N
    c0(j,i) = c0(i,j)
   end do
  end do

 if(present(normalize))then
  if(normalize .eqv. .true.)then
   do i=1,N
    do j=1,N
     c0(i,j) = c0(i,j)/sqrt(c0(i,i)*c0(j,j)) 
    end do
   end do
  end if
 end if
end subroutine

subroutine normalize_matrix(c0)
real,dimension(:,:),pointer :: c0
integer :: i,j,N
   N = size(c0(1,:))
   do i=1,N
    do j=1,N
     c0(i,j) = c0(i,j)/sqrt(c0(i,i)*c0(j,j)) 
    end do
   end do
end subroutine
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine b_factors_ANM(c0,selected_atoms,beta)
use nrutil
real,dimension(:,:),pointer :: c0
type(atom_entry),dimension(:),pointer :: selected_atoms
real,dimension(:),pointer :: beta
!---------------------------------
integer :: i,N
real :: prom,prom1,gamma
real,dimension(1:size(beta)) :: b1
!---------------------------------
 N = size(beta)
 !print *,"Calculating B-factors..."
 do i=1,N !Cycle over super-elements
  beta(i) = c0(i,i)
 end do
 
 print *,"Computing optimal force constant..."
 prom = 0
 do i=1,N
  prom = prom + beta(i)*selected_atoms(i)%beta
  prom1 = prom1 + beta(i)**2
  b1(i) = selected_atoms(i)%beta
  !!
  !print *,beta(i),selected_atoms(i)%beta,prom,prom1
 end do

 gamma = prom/prom1
 beta = beta*8*PI*PI*gamma
 print *,"Force constant = ",gamma," kT"
 print *,"Correlation of B-factors = ",corr(beta,b1)
end subroutine

subroutine make_p_matrix_GNM(c0,p) ! Updated to make computation optimal
real,dimension(:,:),pointer :: c0,p
!--------------------------------
integer :: i,j,N
real :: s1,s2,m1,m2,prod
!--------------------------------
 N = size(c0(1,:))
 print *,"Calculating column-correlations..." 
 do i=1,N
  do j=i,N
   m1 = sum(c0(i,:))/real(N)
   m2 = sum(c0(j,:))/real(N)
   s1 = sum(c0(i,:)**2)/real(N) - m1**2
   s2 = sum(c0(j,:)**2)/real(N) - m2**2
  
   p(i,j) = sum(c0(i,:)*c0(j,:))/real(N) - m1*m2
   p(i,j) = p(i,j) / sqrt(s1*s2)
   
   p(j,i) = p(i,j)
  end do 
 end do
end subroutine

subroutine write_vectors(atoms,filename,vectors,values,num)
type(atom_entry),dimension(:),pointer :: atoms
real,dimension(:),pointer :: values
real,dimension(:,:),pointer :: vectors
character(*) :: filename
integer :: num,i,N

N = size(atoms)
print *,"Writing eigenvectors for ",N," atoms"
open(111,file=trim(filename),action="write")
 do i=1,N
  write(111,"(A2,$)") atoms(i)%chain
 end do
 write(111,*)
 do i=1,N
  write(111,"(I4,$)") atoms(i)%res_num
 end do
 write(111,*)
 do i=3*N-6-num,3*N-6
  !fmt = "(F12.5,"//int_to_str(N*3)//"F12.5)"
  !print *,fmt
  write(111,*) values(i),vectors(:,i)
 end do
close(111)
end subroutine


end module
