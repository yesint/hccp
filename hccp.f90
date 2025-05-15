module hccp
!implicit none
use pdb
use functions
use comm

type prop_type
 real :: interdom_corr 
 real :: dom_corr
 real :: interdom_en
 real :: dom_en
 real :: all_en
 !Domain statistics properties
 real :: seg_num1
 real :: seg_num2
 real :: dom_corr1
 real :: dom_corr2
 real :: dom_size1
 real :: dom_size2
end type

contains

 subroutine do_clustering(c,res,num_natural,interval_natural,num_2,natural_ratio,atoms,do_ISE,printTCL)
 implicit none
 real,dimension(:,:),pointer :: c !Input matrix
 real,dimension(:,:),pointer :: res !Output information
 !character(*),intent(in) :: method
 character(*),optional,intent(inout) :: do_ISE,printTCL
 type(atom_entry),dimension(:),pointer :: atoms
 integer,intent(out) :: num_natural
 real,intent(out) :: natural_ratio
 !---------------------------------------
 integer :: i,j,k,N,old_i,old_num,sz,sz1,sz2,im,jm,cur_proc
 real :: prom,prom1,mer
 real,dimension(0:N_bins) :: bins
 integer,dimension(:),pointer :: cluster,vec2,old_vec
 integer,dimension(:,:),pointer :: cl !Contains info about clusters
 real,dimension(:,:),pointer :: ref_c !Reference matrix
 !---------------------------------------
 integer,dimension(1:size(c(1,:)),1:2) :: segments
 integer :: largest_seg,largest_ind,num_seg,num_clusters,num_big,cur_change
 logical :: bin_empty
 integer :: do_dfire
 real :: mean_internal,mean_external,corr_2,corr_ex
 integer :: interval2,num_2
 
 integer :: interval_natural,interval_cur
 integer,dimension(:),allocatable :: intervals

 integer,dimension(1) :: d1

 integer :: size1,size2,cur
 !---------------------------------------
  res = 0

  if(present(do_ISE) .eqv. .false.)then
   do_ISE="no"
  end if

  if(present(printTCL) .eqv. .false.)then
   printTCL = "no"
  end if

  if(printTCL=="yes")then
   call init_tcl()
  end if

  N = size(c(1,:))
  allocate(cluster(1:N),vec2(1:N),old_vec(1:N))
  allocate(cl(1:N,1:N),ref_c(1:N,1:N))
  allocate(intervals(1:N))
  
  cluster(1:N) = 1 ! N clusters of size 1 are now present
  cl = 0
  do i=1,N
   cl(i,1)=i !Initialize clusters
  end do
  !Find minimum and maximum
  prom  = minval( c, mask=(c/=1) )
  prom1 = maxval( c, mask=(c/=1) )

  !Divide the interval of values into bins
  do i=0,N_bins
   bins(i) = prom1-i*(prom1-prom)/N_bins
  end do
  
 !Start with the first bin
 k = 1
 cur_proc = 1
 do i=1,N
  ref_c(i,1:N) = c(i,1:N)
 end do

 old_i = 0 !Last step was not performed yet
 old_num = 1 !Fake number of clusters
 vec2=0
 old_vec=0
 !call init_tcl()
 !Initialize natural number of clusters count
 interval_natural = 0
 interval_cur = 0
 num_natural = 0
 interval2 = 0
 num_2=0
 intervals = 0
 num_clusters = N
 !Main cycle
 DO WHILE(count(cluster>0)>1 )

  if(mod(k,100)==0) print "(A,I3,A,F10.4,A)",">>> Bin ",k," Value ",bins(k)," ----------------------------------"

  cur_change = 0
  bin_empty = .true. !Assume that bin is empty
  cur = 1 !In order to enter the inner cycle
  do while(cur>0) !Do while everything inside this bin is clustered
   cur = 0 !Nothing clustered yet
   if(old_i/=0)then !If last step was successful then...
    do j=1,N !Update corresponding row and column of c matrix
             !But only if this residues are still in play
     if(cluster(old_i)==1 .and. cluster(j)==1 .and. old_i/=j )then
      call compare_two_mean(old_i,j)
      c(j,old_i) = c(old_i,j)
     end if
    end do
   end if

   !End of matrix update.
   prom = -3  !Fake number
   im = 0
   jm = 0
   !Cycle over all pairs
   out: do i=1,N-1
    do j=i+1,N
     !if we are above current bin and "residues" are in play...
     if(c(i,j)>bins(k) .and. cluster(i)==1 .and. cluster(j)==1 .and. abs(i-j)>0 )then
      !Mark this pair as candidates if their correlation is larger then prev.
      if(c(i,j)>prom)then
       im = i
       jm = j
      end if
     end if
    end do
   end do out
   !If something was found:
   if(im/=0 .and. jm/=0)then
      !Compute the sizes of clusters
      size1 = count(cl(im,:)>0)
      size2 = count(cl(jm,:)>0)
      !Add smaller cluster to larger one
      if(size1>=size2)then
       cl(im,size1+1:size1+size2) = cl(jm,1:size2)
       cluster(jm) = 0 !This cluster is out of play now.
       old_i = im !Cluster i need to be updated next step
      else
       cl(jm,size2+1:size2+size1) = cl(im,1:size1)
       cluster(im) = 0
       old_i = jm !Cluster j need to be updated next step
      end if
      mer = c(im,jm) !Merging was accomplished at mer value of correaltion
      cur = cur + 1 !Increase number of merged pairs
      bin_empty = .false. !Bin is not empty
   end if
  end do !cur
  !We have finished with the current bin
 
  ! Save entry in the stability gap table
  SGtable(k) = count( cluster==1 )

  if(bin_empty .eqv. .false.)then

   !Deal with the natural number of clusters
  
   intervals(num_clusters) = interval_cur
   if(interval_cur >= interval_natural)then
    interval_natural = interval_cur
    num_natural = num_clusters
   end if
   interval_cur = 1

   ! Alternative criterion of the natural number of clusters
   ! span_Ncl - number of largest clusters to compute span (now works only for 2)

   ! span = (cl1_size + cl2_size)/N_bins
   ! span_ratio = cl1_size/cl2_size 
   ! criterion = span * span_ratio    <- this should be maximized
   
   !Mark current clusters
   num_clusters = count( cluster==1 )

   !Fill vec2
   vec2 = 0
   do i=1,N !Cycle over all residues
    if(cluster(i)==1)then !This cluster is present
     sz = count(cl(i,:)>0)
     !Fill vec2
     do j=1,sz
      vec2(cl(i,j))=i
     end do
    end if
   end do !Fill vec2
   
   if(num_clusters>1 .and. trim(do_ISE)=="yes")then
     call ISE()
   end if

   ! Save new entry in the results table
   res(num_clusters,1:N) = vec2(1:N)

  else
   ! The bin is empty - increase interval
   interval_cur = interval_cur + 1
  end if !not empty

  !Print TCL
  if(printTCL=="yes")then
   if(num_clusters/=old_num)then
    cur_proc = cur_proc + 1
    call print_tcl(num_clusters,cur_proc)
    old_num = num_clusters
   end if
  end if

  !Increase current bin
  k = k+1

 END DO
 !See if remaining space sais that there is one cluster
 intervals(1) = N_bins-k
 if(N_bins-k >= interval_natural)then
   interval_natural = N_bins-k
   num_natural = 1
 end if

 d1 = maxloc(intervals)
 print *,"Natural from max:",d1(1)
 print *,"Natural from if :",num_natural

 intervals(d1(1))=0
 d1 = maxloc(intervals)

 print *,"Second natural:",d1(1),intervals(d1(1))
 interval2 = intervals(d1(1))
 num_2 = d1(1)

 natural_ratio = real(interval_natural)/real(interval2)

 do i=1,N
 c(i,1:N) = ref_c(i,1:N) !Restore initial array
 end do

 deallocate(cluster,vec2,ref_c,cl,old_vec,intervals)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 
 CONTAINS

  subroutine sort_it(arr)
  integer,dimension(:),intent(inout) :: arr
  integer :: Nn,i,j,k,prom
  !--------------------------
   Nn = size(arr)
   
   where(arr==0)
    arr=1e6
   end where

   do k=1,Nn-1
   do i=Nn-1,k,-1
    if(arr(i+1)<arr(i) .and. arr(i)>0 .and. arr(i+1)>0)then
     prom = arr(i+1)
     arr(i+1) = arr(i)
     arr(i) = prom
    end if
   end do
   end do

   where(arr==1e6)
    arr=0
   end where

  end subroutine

  subroutine sort_out_0(arr)
  integer,dimension(:),intent(inout) :: arr
  integer :: Nn,i,j,k,prom
  !--------------------------
   Nn = size(arr)
   do k=Nn,2,-1
   do i=2,k
    if(arr(i-1)<arr(i) )then
     prom = arr(i-1)
     arr(i-1) = arr(i)
     arr(i) = prom
    end if
   end do
   end do
  end subroutine

 
  subroutine compare_two_mean(i,j)
  integer,intent(in) :: i,j
  integer :: sz1,sz2,k,l
  real :: prom,prom1,prom2,prom_c,cm1,cm2
  !------------------------
  !Determine size of the compared clusters
    sz1 = count(cl(i,:)>0)
    sz2 = count(cl(j,:)>0)
    !print "(<N>I4)",cl(i,:)
    !print "(<N>I4)",cl(j,:)
    prom = 0
    do k=1,sz1
     do l=1,sz2
      !print *,"$",i,j,k,l
      prom = prom + ref_c(cl(i,k),cl(j,l))
     end do
    end do
    prom = prom/real(sz1*sz2)
    c(i,j) = prom
  end subroutine

!-----------------------------------------------------------------------------------

  subroutine verify_all()
  integer :: i,j,k,sz
  integer,dimension(1:N) :: h
   print *,"VERIFYING..."
   h=0
   do i=1,N
    if(cluster(i)/=0)then !If cluster is in play
     sz = count(cl(i,:)>0)
     do j=1,sz
      h(cl(i,j)) = h(cl(i,j))+1
     end do
    end if
   end do

   if(count(h/=1)>0)then
    print *,"<ERROR>"
    do i=1,N
     print "(I4,$)",i
    end do
    print *
    fmt = "("//int_to_str(N)//"I4)"
    print fmt,h
   end if
  end subroutine

!/------------------------------------------------------------------------------------
  subroutine ISE() !Intercalating segments removing
  integer,dimension(1:3) :: left, middle, right
  integer,dimension(:),pointer :: cl_IS,cl_rest,cl_out
  integer :: i,j,k,cur_ind,cur_val,sz_IS,sz_all,sz_out
  real :: c1,c2
  !----------------------------
    integer,dimension(1:N) :: tt
  !----------------------------

   !!
   !print *,"==== ISE ==== old_i = ",old_i
   !!

   allocate(cl_IS(1:N),cl_rest(1:N))
   !Find left
78 cur_ind = 1
   left =0
   middle=0
   right=0
   cl_IS=0
   cl_rest=0
   cur_val = vec2(cur_ind)
   left(1) = cur_ind
   left(3) = cur_val
   do while(vec2(cur_ind)==cur_val .and. cur_ind<=N) 
    cur_ind = cur_ind + 1
   end do
   left(2) = cur_ind-1
   if(left(2)==N) goto 511 !The segment spans the whole protein
   !print *,"left",left(1:3)
   !Find middle
   middle(1) = cur_ind
   cur_val = vec2(cur_ind)
   middle(3) = cur_val
   do while(vec2(cur_ind)==cur_val .and. cur_ind<=N) 
    cur_ind = cur_ind + 1
   end do
   middle(2) = cur_ind-1
   if(middle(2)==N) goto 511 !The segment spans the whole protein
   !print *,"middle",middle(1:3)
   
   !Now enter the cycle
   DO
    !Finding right
    right(1) = cur_ind
    cur_val = vec2(cur_ind)
    right(3) = cur_val
    do while(vec2(cur_ind)==cur_val .and. cur_ind<=N) 
     cur_ind = cur_ind + 1
    end do
    right(2) = cur_ind-1
    !All three segments are found
    if(right(3)==left(3))then
     !Candidate for IS is found
     sz_all = count(cl(middle(3),:)>0)
     sz_IS = middle(2) - middle(1) + 1
     if(sz_IS < sz_all .and. sz_IS<=10)then
      !This is really IS
      
      !!
      !print *,"**** IS FOUND ****"
      !do i=1,N
      ! tt(i) = i
      !end do
      !print *,"!-- current vec2:"
      !print "(<N>I4)",tt
      !print "(<N>I4)",vec2
      !print *,"IS found:",middle(1),middle(2)," value",middle(3)," size",sz_IS," all",sz_all
      !print *,"cl, cl_IS, cl_rest, cl_out:"
      !print "(<N>I4)",cl(middle(3),1:)
      !!

      !Fill cl_is
      do i=1,sz_IS
       cl_IS(i) = middle(1)+i-1
      end do
      !Fill cl_rest. Include all except those from cl_IS
      k=0
      do i=1,sz_all
       if(cl(middle(3),i)<middle(1) .or. cl(middle(3),i)>middle(2) )then
        k = k+1
        cl_rest(k) = cl(middle(3),i)
       end if
      end do
      !Make cl_out
      cl_out=>cl(left(3),:)
      sz_out = count(cl_out>0)

      !!
      ! print "(<N>I4)",cl_IS(1:)
      ! print "(<N>I4)",cl_rest(1:)
      ! print "(<N>I4)",cl_out(1:)
      !!

      !Calculate correlations
      c1 = compare_2cl(cl_IS,cl_rest,ref_c)
      c2 = compare_2cl(cl_IS,cl_out,ref_c)

      !!
      ! print *,"Corr with same :",c1
      ! print *,"Corr with outer:",c2
      !!

      !If correlation with the outer cluster is stronger - update clusters
      if(c2>c1)then
       !!
       !  print *,"Update initiated:"
       !  call verify_all()       !1
       !!

       !Update vec2
       vec2(middle(1):middle(2)) = left(3)
       !Merge cl_IS with cl_out
       cl_out(sz_out+1:sz_out+sz_IS) = cl_IS(1:sz_IS)

       !print *, "Out:"
       !print "(<N>I4)",cl_out(1:N)

       !Delete cl_IS from middle cluster
       do i=1,sz_all
        if(cl(middle(3),i)>=middle(1) .and. cl(middle(3),i)<=middle(2) )then
         cl(middle(3),i) = 0
        end if
       end do

       !print *,"With Zeroes:"
       !print "(<sz_all>I4)",cl(middle(3),1:sz_all)
       !Sort out zeros
       call sort_it(cl(middle(3),:))

       !!
       ! print "(<N>I4)",cl(middle(3),1:)
       ! print "(<N>I4)",cl_out(1:)
       !!

       !We changed two clusters, so update correlation matrix
       do j=1,N !Update cluster cl_IS:
        call compare_two_mean(middle(3),j)
        c(j,middle(3)) = c(middle(3),j)
        call compare_two_mean(left(3),j)
        c(j,left(3)) = c(left(3),j)
       end do

       print *,"INFO> Intercalating segment corrected."
       !call verify_all()         !6
       !print *,"Starting over"
       goto 78

      end if

     end if
    end if
    !End of IS calculations
    left = middle
    middle = right

    if(right(2)==N) exit !Exit the cycle if last segment is found
   END DO
   deallocate(cl_IS,cl_rest)

   !Update vec2
   !Mark current clusters
   num_clusters = count( cluster==1 )

   !Fill vec2
   vec2 = 0
   do i=1,N !Cycle over all residues
    if(cluster(i)==1)then !This cluster is present
     sz = count(cl(i,:)>0)
     !Fill vec2
     do j=1,sz
      vec2(cl(i,j))=i
     end do
    end if
   end do !Fill vec2
      
511 end subroutine

!/------------------------------------------------------------------------------------

!=======================================================
subroutine init_tcl()
!Initialize tcl script
open(111,file="color.tcl",action="write")
 write(111,*) "set num_proc 0"
 write(111,*) "set cur 0"
 write(111,*) "proc next {} {"
 write(111,*) " global cur"
 write(111,*) " incr cur 1"
 write(111,*) " Show$cur"
 write(111,*) "}"
 write(111,*) "proc prev {} {"
 write(111,*) " global cur"
 write(111,*) " set cur [expr $cur-1]"
 write(111,*) " Show$cur"
 write(111,*) "}"
 write(111,*) "proc last {} {"
 write(111,*) " global num_proc"
 write(111,*) " Show$num_proc"
 write(111,*) "}"
 write(111,*) "proc goto {num} {"
 write(111,*) " global cur"
 write(111,*) " set cur $num"
 write(111,*) " Show$cur"
 write(111,*) "}"
close(111)
end subroutine


subroutine print_tcl(m,step) !Print tcl script for given Nber of clusters m
integer,intent(in) :: m,step 
integer,dimension(1:m) :: clusters
integer :: pr,color,i,j,k,res
character(3) :: ch
 !Update the vec in order to mark clusters sequentially 1,2,3...
 clusters = 0
 clusters(1) = vec2(1)
 do i=1,N
  k=0
  do j=1,m
   if(vec2(i)/=clusters(j) .and. clusters(j)/=0)then
    k = k+1
    pr = vec2(i)
   end if
   if(vec2(i)==clusters(j))then
    k = 0
    exit
   end if
  end do
  if(k>0)then
   clusters(k+1)=pr
  end if
 end do
 do i=1,m
  where(vec2==clusters(i))
   vec2=i
  end where
 end do

 !vec2 = mark_sequencially(vec2,m)

 !Now form tcl script
 write(ch,"(I3)") step !Convert step to string
 open(111,file="color.tcl",action="write",access="append")
 write(111,*) "#-- Number of clusters ",m,"  --------------------------------"
 write(111,*) "incr num_proc 1"
 write(111,*) "proc Show",adjustl(ch)," {} {"
 write(111,*) "global cur num_proc"
 write(111,*) "set numRep [molinfo top get numreps]"
 write(111,*) "incr numRep 1"
 write(111,*) "for {set i 0} {$i < $numRep} {incr i 1} {"
 write(111,*) " mol delrep 0 top"
 write(111,*) "}"
 write(111,*) "mol rep NewCartoon"
 write(111,*) "mol color colorid 0"
 write(111,*) 'mol selection "backbone"'
 write(111,*) "mol addrep top"
 color = 1
 do i=1,m
  if(count(vec2==i)>5)then
  write(111,*) "mol color colorid ",color
  color = color + 1
  write(111,"(a,$)") 'mol selection "backbone and ('
  do j=1,N
   if(vec2(j)==i)then
    if(atoms(j)%chain/=" ")then
     write(111,"(a,i5,a,a,a,$)") "(resid=",atoms(j)%res_num," and chain ",atoms(j)%chain," ) or "
    else
     write(111,"(a,i5,a,$)") "(resid=",atoms(j)%res_num," ) or "
    end if
   end if
  end do
  write(111,"(a)") '0=1)"'
  write(111,*) "mol addrep top"
  end if
 end do
 write(111,*) 'puts "$cur of $num_proc"'
 !Print cluster as a range of residues
 pr=vec2(1)
 write(111,"(a,I3,a,a,I3,a,$)") "# (",pr,"|",atoms(1)%chain,atoms(1)%res_num,":"
 do i=1,N
  if(vec2(i)/=pr)then
   pr=vec2(i)
   write(111,"(a,I3,a)") atoms(i-1)%chain,atoms(i-1)%res_num,")"
   write(111,"(a,I3,a,a,I3,a,$)") "# (",pr,"|",atoms(i)%chain,atoms(i)%res_num,":"
  end if
 end do
 write(111,"(a,I3,a)") atoms(N)%chain,atoms(N)%res_num,")"
 write(111,*) 

 write(111,*) "}"
 close(111)

end subroutine

end subroutine !hccp

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

 subroutine corr_and_energy(vec,atoms,ref_c,n_dom,prop)
 integer,dimension(:),pointer :: vec
 type(atom_entry),dimension(:),pointer :: atoms,pot_atoms
 type(prop_type) :: prop
 integer :: n_dom
 real,dimension(:,:),pointer :: ref_c
 !----------------------
 real,dimension(:,:,:),pointer :: pot
 integer,dimension(:,:),allocatable,target :: cl
 integer :: i,j,k,N,num_seg,sz
 integer,dimension(:),pointer :: cl1,cl2
 real :: prom,prom2,pr,prom3
 integer,dimension(1:25,1:2) :: segments

 !integer,dimension(1:N) :: hinges
 !integer :: num_hinges
 !----------------------
  call select_atoms(atoms,pot_atoms,"SCG",N)
  allocate(cl(1:n_dom,1:N))
  !Load dfire potential
  print *,"Reading DFIRE potential "
  allocate(pot(1:20,1:20,1:20))
  call read_dfire("dfire_scm.pot",pot)
  !Parse domain each to separate cl
  cl = 0
  do i=1,n_dom
   k=0
   do j=1,N
    if(vec(j)==i)then
     k = k+1
     cl(i,k)=j
    end if
   end do
  end do

  !Try to find hinges
  !To do so find points, which are "anchors" for possible hinges
  

  if(n_dom>1)then
   if(printDB=="yes")then
   open(111,file="assignment.dat",action="write",access="append")
   do i=1,n_dom
    !Begin extracting segments for current cluster
    sz = count(cl(i,:)>0)
    write(111,"(A,I4,A,I4)") "  Domain",i," Size: ",sz
    segments = 0
    num_seg = 1
    segments(1,1)=cl(i,1)
    !call sort_it(cl(i,:))
    do j=1,sz-1 !Cycle over all residues in cluster
     if( cl(i,j+1)-cl(i,j)/=1 )then
      segments(num_seg,2) = cl(i,j)
      num_seg = num_seg+1
      segments(num_seg,1) = cl(i,j+1)
     end if
    end do
    segments(num_seg,2)=cl(i,sz)
    !Print info about segments
    do j=1,num_seg
     write(111,"(A,I3,A,I3,A,$)") "  [",segments(j,1),":",segments(j,2),"] "
     !print "(A,I3,A,I3,A,$)","  [",segments(j,1),":",segments(j,2),"] "
    end do

    !Fill some statistics (only meaningfull for forced 2 domains)
    if(i==1)then
      prop%seg_num1 = num_seg
      prop%dom_size1 = real(sz)/real(N)
    end if
    if(i==2)then
      prop%seg_num2 = num_seg
      prop%dom_size2 = real(sz)/real(N)
    end if

    write(111,*) ""
    !print *,""
   end do
   close(111)
   end if !printDB
  end if

  !!
   ! print "(<N>I4)",vec(1:N)
   ! do i=1,n_dom
   !  print "(<N>I4)",cl(i,1:N)
   ! end do
  !!
  
  !print *,"Calculate domain correlations..."
  prom = 0
  do i=1,n_dom
   cl1 => cl(i,1:N)
   prom2 = compare_1cl(cl1,ref_c)
   prom = prom + prom2

   !Fill some more statistics
   if(i==1)then
      prop%dom_corr1 = prom2
   end if
   if(i==2)then
      prop%dom_corr2 = prom2
   end if

  end do
  prop%dom_corr = prom/real(n_dom)

  !Calculate interdomain correlations (cycle over domain pairs)
  if(n_dom>1)then
   prom = 0
   do i=1,n_dom-1
    do j=i+1,n_dom
     cl1 => cl(i,1:N)
     cl2 => cl(j,1:N)
     prom2 = compare_2cl(cl1,cl2,ref_c)
     !print *,"Pair ",i,j,prom2
     prom = prom + prom2
    end do
   end do
   prop%interdom_corr = prom/real((n_dom**2-n_dom)/2)
  else
   prop%interdom_corr = 0.0
  end if

  !Calculate domain and interdomain energy
  prom = dfire_energy(pot_atoms,pot)
  prop%all_en = prom/real(N)
  prom2 = 0
  prom3 = 0
  do i=1,n_dom
   cl1 => cl(i,1:N)
   !print "(<N>I4)",cl1(1:N)
   pr = dfire_energy_index(pot_atoms,cl1,pot)
   prom2 = prom2 + pr
   prom3 = prom3 + pr/count(cl1>0)
  end do
  prop%interdom_en = (prom-prom2)/real(N)
  prop%dom_en = prom3/n_dom

  deallocate(pot,cl,pot_atoms)
 end subroutine

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  function compare_2cl(cl1,cl2,ref_c)
  integer,dimension(:),pointer :: cl1,cl2
  real,dimension(:,:),pointer :: ref_c
  real :: compare_2cl
  integer :: sz1,sz2,k,l
  real :: prom
  !------------------------
  !Determine size of the compared clusters
    sz1 = count(cl1(:)>0)
    sz2 = count(cl2(:)>0)
    !print "(<N>I4)",cl(i,:)
    !print "(<N>I4)",cl(j,:)
    prom = 0
    do k=1,sz1
     do l=1,sz2
      prom = prom + ref_c(cl1(k),cl2(l))
     end do
    end do
    prom = prom/real(sz1*sz2)
    compare_2cl = prom
  end function

  function compare_1cl(cl1,ref_c)
  integer,dimension(:),pointer :: cl1
  real,dimension(:,:),pointer :: ref_c
  real :: compare_1cl
  integer :: sz1,k,l
  real :: prom
  !------------------------
  !Determine size of the compared clusters
    sz1 = count(cl1(:)>0)
    prom = 0
    do k=1,sz1-1
     do l=k+1,sz1
      prom = prom + ref_c(cl1(k),cl1(l))
     end do
    end do
    prom = prom/real((sz1**2-sz1)/2)
    compare_1cl = prom
  end function

function mark_sequencially(vec2,m) !m - number of clusters
integer,dimension(:),pointer :: vec2
integer,dimension(1:size(vec2)) :: mark_sequencially
integer :: i,k,j,N,m,pr
integer,dimension(1:m) :: clusters
!--------------------------
 !print *,"MARK_SEQUENCIALLY enterd. m=",m
 N = size(vec2)
 !print "(<N>I4)", vec2
 clusters = 0
 clusters(1) = vec2(1)
 do i=1,N
  k=0
  do j=1,m
   if(vec2(i)/=clusters(j) .and. clusters(j)/=0)then
    k = k+1
    pr = vec2(i)
   end if
   if(vec2(i)==clusters(j))then
    k = 0
    exit
   end if
  end do
  if(k>0)then
   clusters(k+1)=pr
  end if
 end do
 !Re-mark clusters sequencially
 do i=1,m
  where(vec2==clusters(i))
   mark_sequencially=i
  end where
 end do
 !print "(<m>I4)", clusters
 !print "(<N>I4)", mark_sequencially
end function

end module
