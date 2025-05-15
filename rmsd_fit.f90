module rmsd_fit
!Performs RMSD fit of two structures
!The code is hacked from GROMACS 3.1.4 from the file do_fit.c
!Do not ask me how it works - I've no idea about algorithm.
!Translated from C and thus looks ugly :( but works perfectly :))
use pdb
use comm
use eigen

contains
!Brings atoms' center of masses to 0
subroutine center(atoms)
type(atom_entry),dimension(:),pointer :: atoms
!-----------------------
integer :: i,N
real,dimension(1:3) :: cm
real :: total_mass
!-----------------------
 N = size(atoms)
 cm = 0
 total_mass = 0
 do i=1,N
  if(atoms(i)%mass/=0)then
   cm(1:3) = atoms(i)%xyz(1:3) * atoms(i)%mass
   total_mass = total_mass + atoms(i)%mass
  else
   cm(1:3) = atoms(i)%xyz(1:3) 
   total_mass = total_mass + 1
  end if
 end do
 cm = cm/total_mass
 
 do i=1,N
  atoms(i)%xyz(1:3) = atoms(i)%xyz(1:3) - cm(1:3)
 end do
end subroutine

!Calculates rotation matrix
subroutine rot_matr(ref_atoms,atoms,R)
type(atom_entry),dimension(:),pointer :: atoms,ref_atoms
real,dimension(1:3,1:3) :: R
!---------------------------
integer :: i,j,N,j1,j2
real,dimension(1:6,1:6) :: omega,om
real,dimension(1:3,1:3) :: u,vh,vk
real,dimension(1:6) :: d
real :: xnr,xpc,max_d
integer :: nrot,index
!---------------------------
 N = size(atoms)
 d = 0
 omega = 0
 om = 0
 !Calculate the matrix U
 u = 0
 do i=1,N
  do j1=1,3
   xpc = ref_atoms(i)%xyz(j1)
   do j2=1,3
    xnr = atoms(i)%xyz(j2)
    u(j1,j2) = u(j1,j2) + xnr*xpc*atoms(i)%mass
   end do
  end do
 end do

 !Construct omega
 do i=1,6
  do j=1,i
   if(i>3 .and. j<=3)then
    omega(i,j) = u(i-3,j)
    omega(j,i) = omega(i,j)
   else
    omega(i,j) = 0
    omega(j,i) = 0
   end if
  end do
 end do

 !Finding eigenvalues of omega
 call jacobi(omega,d,om,nrot)
 !print *,"Digonalized with ",nrot," rotations"

 do j=1,2
    max_d=-1000;
    do i=1,6
      if(d(i)>max_d)then
        max_d=d(i)
        index=i
      end if
    end do
    d(index)=-10000;
    do i=1,3
      vh(j,i)=sqrt(2.0)*om(i,index);
      vk(j,i)=sqrt(2.0)*om(i+3,index);
    end do
 end do

 ! Calculate the last eigenvector as the outer-product of the first two.
 ! This insures that the conformation is not mirrored and
 ! prevents problems with completely flat reference structures.

  vh(3,:) = oprod(vh(1,:),vh(2,:)) 
  vk(3,:) = oprod(vk(1,:),vk(2,:))

 !Determine R matrix
 do i=1,3
  do j=1,3
   R(i,j) = vk(1,i)*vh(1,j) + vk(2,i)*vh(2,j) + vk(3,i)*vh(3,j);
  end do
 end do

end subroutine

function oprod(a,b)
real,dimension(1:3) :: a,b,oprod
  oprod(1)=a(2)*b(3)-a(3)*b(2);
  oprod(2)=a(3)*b(1)-a(1)*b(3);
  oprod(3)=a(1)*b(2)-a(2)*b(1);
end function

!Subroutine to actually do fitting
! Fits atoms to ref_atoms in the same way as VMD does
! atoms and ref_atoms should be of the same size
! The coordinates of both atoms and ref_atoms are ALTERED AND NOT RESTORED!!!
subroutine do_fit(ref_atoms,atoms)
type(atom_entry),dimension(:),pointer :: atoms,ref_atoms
!------------------------------
integer :: i,j,N,m,r,c
real,dimension(1:3,1:3) :: RR
real,dimension(1:3) :: old
!------------------------------
 N = size(atoms)
 call center(ref_atoms)
 call center(atoms)
 call rot_matr(ref_atoms,atoms,RR)
 !Apply rotation
 do j=1,N
    do m=1,3
     old(m)=atoms(j)%xyz(m)
    end do
    do r=1,3
      atoms(j)%xyz(r)=0
      do c=1,3
        atoms(j)%xyz(r) = atoms(j)%xyz(r) + RR(r,c)*old(c)
      end do
    end do
  end do
end subroutine

end module