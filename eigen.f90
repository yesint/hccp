module eigen

contains

SUBROUTINE jacobi(a,d,v,nrot)
USE nrtype; USE nrutil, ONLY : assert_eq,get_diag,nrerror,unit_matrix,upper_triangle
IMPLICIT NONE
INTEGER(I4B), INTENT(OUT) :: nrot
REAL(SP), DIMENSION(:), INTENT(OUT) :: d
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
!Computes all eigenvalues and eigenvectors of a real symmetric N .N matrix a. On output,
!elements of a above the diagonal are destroyed. d is a vector of length N that returns the
!eigenvalues of a. v is an N . N matrix whose columns contain, on output, the normalized
!eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
INTEGER(I4B) :: i,ip,iq,n,ii,jj
REAL(SP) :: c,g,h,s,sm,t,tau,theta,tresh
REAL(SP), DIMENSION(size(d)) :: b,z
 n=assert_eq((/size(a,1),size(a,2),size(d),size(v,1),size(v,2)/),'jacobi')
 call unit_matrix(v(:,:)) !Initialize v to the identity matrix.
 b(:)=get_diag(a(:,:)) !Initialize b and d to the diagonal of a. 
 d(:)=b(:)
 z(:)=0.0 !This vector will accumulate terms of the form tapq as in eq. (11.1.14). 
 nrot=0
 do i=1,50
  !sm=sum(abs(a),mask=upper_triangle(n,n)) !Sum o.-diagonal elements.
  sm = 0
  do ii=1,size(d)
   do jj=ii+1,size(d)
    sm = sm + abs(a(ii,jj))
   end do
  end do
  if (sm == 0.0) RETURN
  !  The normal return, which relies on quadratic convergence to machine under.ow.
  tresh=merge(0.2_sp*sm/n**2,0.0_sp, i < 4 )
  !On the .rst three sweeps, we will rotate only if tresh exceeded.
  do ip=1,n-1
   do iq=ip+1,n
    g=100.0_sp*abs(a(ip,iq))
    !After four sweeps, skip the rotation if the o.-diagonal element is small.
    if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
    .and. (abs(d(iq))+g == abs(d(iq)))) then
      a(ip,iq)=0.0
    else if (abs(a(ip,iq)) > tresh) then
      h=d(iq)-d(ip)
      if (abs(h)+g == abs(h)) then
       t=a(ip,iq)/h !t = 1/(2.)
      else
       theta=0.5_sp*h/a(ip,iq) !Equation (11.1.10).
       t=1.0_sp/(abs(theta)+sqrt(1.0_sp+theta**2))
       if (theta < 0.0) t=-t
      end if
      c=1.0_sp/sqrt(1+t**2)
      s=t*c
      tau=s/(1.0_sp+c)
      h=t*a(ip,iq)
      z(ip)=z(ip)-h
      z(iq)=z(iq)+h
      d(ip)=d(ip)-h
      d(iq)=d(iq)+h
      a(ip,iq)=0.0
      call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
      !Case of rotations 1 öœ j < p.
      call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
      !Case of rotations p < j < q.
      call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
      !Case of rotations q < j öœ n.
      call jrotate(v(:,ip),v(:,iq))
      nrot=nrot+1
     end if
    end do
   end do
   b(:)=b(:)+z(:)
   d(:)=b(:) !Update d with the sum of tapq,
   z(:)=0.0 !and reinitialize z.
  end do
  call nrerror("too many iterations in jacobi")

contains

 SUBROUTINE jrotate(a1,a2)
 REAL(SP), DIMENSION(:), INTENT(INOUT) :: a1,a2
 REAL(SP), DIMENSION(size(a1)) :: wk1
  wk1(:)=a1(:)
  a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
  a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
 END SUBROUTINE jrotate

END SUBROUTINE jacobi

SUBROUTINE eigsrt(d,v)
USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,swap
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: v
!Given the eigenvalues d and eigenvectors v as output from jacobi (11.1) or tqli (11.3),
!this routine sorts the eigenvalues into descending order, and rearranges the columns of v
!correspondingly. The method is straight insertion.
INTEGER(I4B) :: i,j,n
n=assert_eq(size(d),size(v,1),size(v,2),'eigsrt')
do i=1,n-1
 j=imaxloc(d(i:n))+i-1
 if (j /= i) then
  call swap(d(i),d(j))
  call swap(v(:,i),v(:,j))
 end if
end do
END SUBROUTINE eigsrt
!=====================================================

SUBROUTINE tred2(a,d,e,novectors)
USE nrtype; USE nrutil, ONLY : assert_eq,outerprod
IMPLICIT NONE
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
INTEGER(I4B) :: i,j,l,n,ii,jj
REAL(SP) :: f,g,h,hh,scale
REAL(SP), DIMENSION(size(a,1)) :: gg
LOGICAL(LGT) :: yesvec
 n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2')
 if (present(novectors)) then
  yesvec=.not. novectors
 else
  yesvec=.true.
 end if
 do i=n,2,-1
  l=i-1
  h=0.0
  if (l > 1) then
   scale=sum(abs(a(i,1:l)))
   if (scale == 0.0) then 
    e(i)=a(i,l)
   else
    a(i,1:l)=a(i,1:l)/scale 
    h=sum(a(i,1:l)**2)
    f=a(i,l)
    g=-sign(sqrt(h),f)
    e(i)=scale*g
    h=h-f*g 
    a(i,l)=f-g 
    if (yesvec) a(1:l,i)=a(i,1:l)/h 
     do j=1,l 
      e(j)=(dot_product(a(j,1:j),a(i,1:j)) &
      +dot_product(a(j+1:l,j),a(i,j+1:l)))/h
     end do
     f=dot_product(e(1:l),a(i,1:l))
     hh=f/(h+h) 
     e(1:l)=e(1:l)-hh*a(i,1:l)
     do j=1,l 
      a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
     end do
    end if
   else
    e(i)=a(i,l)
   end if
   d(i)=h
  end do
  if (yesvec) d(1)=0.0
  e(1)=0.0
  do i=1,n 
   if (yesvec) then
    l=i-1
    if (d(i) /= 0.0) then
     gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
     !!!
!     print *,"Trying outerprod ",i,"..."
     !!!
     !a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
     	 do ii=1,l
	  do jj=1,l
           a(ii,jj)=a(ii,jj)-a(ii,i)*gg(jj)
	   !outerprod_r(i,j) = a(i)*b(j)
	  end do
	 end do
     !!!
!     print *,"Done outerprod."
     !!!
    end if
    d(i)=a(i,i)
    a(i,i)=1.0 
    a(i,1:l)=0.0
    a(1:l,i)=0.0
   else
    d(i)=a(i,i)
   end if
  end do
END SUBROUTINE tred2
!=======================================================
SUBROUTINE tqli(d,e,z)
USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
USE nr, ONLY : pythag
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
INTEGER(I4B) :: i,iter,l,m,n,ndum
REAL(SP) :: b,c,dd,f,g,p,r,s
REAL(SP), DIMENSION(size(e)) :: ff
 n=assert_eq(size(d),size(e),'tqli: n')
 if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
  e(:)=eoshift(e(:),1) 
  do l=1,n
   iter=0
   iterate: do
   do m=l,n-1 
    dd=abs(d(m))+abs(d(m+1))
    if (abs(e(m))+dd == dd) exit
   end do
   if (m == l) exit iterate
   if (iter == 30) call nrerror('too many iterations in tqli')
    iter=iter+1
    g=(d(l+1)-d(l))/(2.0_sp*e(l)) 
    r=pythag1(g,1.0_sp)
    g=d(m)-d(l)+e(l)/(g+sign(r,g)) 
    s=1.0
    c=1.0
    p=0.0
    do i=m-1,l,-1 
     f=s*e(i)
     b=c*e(i)
     r=pythag1(f,g)
     e(i+1)=r
     if (r == 0.0) then 
      d(i+1)=d(i+1)-p
      e(m)=0.0
      cycle iterate
     end if
     s=f/r
     c=g/r
     g=d(i+1)-p
     r=(d(i)-g)*s+2.0_sp*c*b
     p=s*r
     d(i+1)=g+p
     g=c*r-b
     if (present(z)) then
      ff(1:n)=z(1:n,i+1)
      z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
      z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
     end if
    end do
    d(l)=d(l)-p
    e(l)=g
    e(m)=0.0
   end do iterate
  end do
END SUBROUTINE tqli

FUNCTION pythag1(a,b)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: a,b
REAL(SP) :: pythag1
!Computes (a2 + b2)1/2 without destructive under.ow or over.ow.
REAL(SP) :: absa,absb
absa=abs(a)
absb=abs(b)
if (absa > absb) then
pythag1=absa*sqrt(1.0_sp+(absb/absa)**2)
else
if (absb == 0.0) then
pythag1=0.0
else
pythag1=absb*sqrt(1.0_sp+(absa/absb)**2)
end if
end if
END FUNCTION pythag1


end module