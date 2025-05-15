program hccp_classic
 
 use pdb
 use functions
 use nrtype
 use hccp
 use comm
 use input_parser

 implicit none

 type(atom_entry),dimension(:),pointer :: atoms,selected_atoms
 real,dimension(:),pointer :: values
 real,dimension(:,:),pointer :: vectors
 integer :: num_atoms,N
 real :: rc
 real,dimension(:,:,:),pointer :: pot
 character(3),dimension(1:20) :: hash
 character(200) :: sel_text,potname
 !-------------------------------------
 integer :: x,y,x0,y0
 integer,dimension(1) :: dum1
 integer :: cur,i,j,k,index,dir,k1,old_dir,cur_dir
 real,dimension(:,:,:),allocatable :: conf
 real :: scale
 real,dimension(:,:),pointer :: ref
 integer,dimension(:),pointer :: vec
 real,dimension(1:100) :: energies
 !-------------------------------------
 real :: en,rms,old_rms,gamma,prom,prom1
 integer :: active
 character(80) :: ch,method
 character(200) :: matrix_filename
 character(1) :: chain
 
 real,dimension(:,:),pointer :: c0,p
 real,dimension(:),pointer :: beta,b1
 type(prop_type) :: properties
 !-------------------------------------
 print *,"============================================="
 print *,"            HCCP Finder, V.1.5"
 print *,"     (C) Semen Yesylevskyy, 2005-2006"
 print *,""
 print *," Please visit "
 print *," http://www.geocities.com/yesint3/hccp.html"
 print *," for information about usage, citations etc."
 print *,"============================================="

 call parse_file("<STDIN>")
 
 call define_keyword(name="file",group="pdb",data_type="string",variable=filename)
 call define_keyword(name="chain",group="pdb",data_type="string",variable=chain)

 call define_keyword(name="method",group="gnm",data_type="string",variable=method,options=[character(len=11) :: "GNM","ANM","read_matrix"])
 call define_keyword(name="selection",group="gnm",data_type="string",variable=sel_text)
 call define_keyword(name="cut_off",group="gnm",data_type="real",variable=rc)
 call define_keyword(name="matrix_file",group="gnm",data_type="string",variable=matrix_filename,default_value="")

 call define_keyword(name="ISE",group="hccp",data_type="string",variable=doISE,options=[character(len=3) :: "yes","no"])
 call define_keyword(name="force_two_domains",group="hccp",data_type="string",variable=force2domains,options=[character(len=3) :: "yes","no"])

 call define_keyword(name="print_tcl",group="output",data_type="string",variable=printTCL,options=[character(len=3) :: "yes","no"])
 call define_keyword(name="print_code_string",group="output",data_type="string",variable=printVEC,options=[character(len=3) :: "yes","no"])
 call define_keyword(name="print_c_as_line",group="output",data_type="string",variable=printCline,options=[character(len=3) :: "yes","no"])
 call define_keyword(name="print_p_as_line",group="output",data_type="string",variable=printPline,options=[character(len=3) :: "yes","no"])
 call define_keyword(name="print_g_as_line",group="output",data_type="string",variable=printGline,options=[character(len=3) :: "yes","no"])
 call define_keyword(name="print_db",group="output",data_type="string",variable=printDB,options=[character(len=3) :: "yes","no"])
 call define_keyword(name="print_mobility",group="output",data_type="string",variable=printMobility,options=[character(len=3) :: "yes","no"])
 call define_keyword(name="eigen_for_mob",group="output",data_type="integer",variable=numVecs)

 call validate_keywords()

 !------------------------------------------------------------------------

 call read_pdb(trim(filename),atoms,num_atoms,chain)
 call select_atoms(atoms,selected_atoms,trim(sel_text),N)

 if(trim(method)=="GNM" .or. trim(method)=="read_matrix")then
  allocate(vectors(1:N,1:N))
  allocate(values(1:N))
 else !ANM
  allocate(vectors(1:3*N,1:3*N))
  allocate(values(1:3*N))
 end if

 allocate(c0(1:N,1:N), p(1:N,1:N))
 allocate(beta(1:N),b1(1:N))

 if(trim(method)=="read_matrix")then
  print *," Reading external correlation matrix from ",trim(matrix_filename),"..."
  open(111,file=matrix_filename,action="read")
   do i=1,N
    read(111,*) c0(i,1:N)
    !print "(<N>F24.5)",c0(i,1:N)
   end do
  close(111)
  goto 876
 end if

 if(trim(method)=="GNM")then
  call do_gnm(N,rc,selected_atoms,values,vectors)
 else
  call do_anm(N,rc,selected_atoms,values,vectors)
  call write_vectors(selected_atoms,"eigen_VMD.dat",vectors,values,20)
  open(111,file="first.dat",action="write")
  do i=1,N
   write(111,"(3F12.5)") vectors(3*(i-1)+1,3*N-6), vectors(3*(i-1)+2,3*N-6), vectors(3*(i-1)+3,3*N-6)
  end do
  close(111)
 end if

 !Print mobilities if needed
 

 if(printMobility=="yes")then
  call make_corr_matrix_GNM(vectors,values,c0,.false.)
  call b_factors_ANM(c0,selected_atoms,beta)
  !sqr of the eigenvectors added together
  b1 = 0
  do i=N-1-numVecs,N-1
   b1(1:N) = b1(1:N) + vectors(:,i)**2/values(i)
  end do
  !Actually print the mobilities
  open(111,file="Mobilities.dat",action="write")
   write(111,*) "resNum resName beta mob_all mob"
   do i=1,N
    write(111,"(I7,A5,3F12.5)") selected_atoms(i)%res_num, selected_atoms(i)%res_name, selected_atoms(i)%beta, beta(i), b1(i)
   end do
  close(111)
 end if

 if(trim(method)=="GNM")then
  call make_corr_matrix_GNM(vectors,values,c0,.true.)
 else
  call make_corr_matrix_ANM(vectors,values,c0,.true.)
 end if

 !876 if(trim(method)=="read_matrix") call normalize_matrix(c0)
 876 call make_p_matrix_GNM(c0,p)

 if(printCline=="yes")then
  open(111,file="c_string.dat",action="write",access="append")
   fmt = "(A,"//int_to_str(N*N)//"F12.5)"
   write(111,fmt) trim(filename)//" ",c0
  close(111)
 end if

 if(printPline=="yes")then
  open(111,file="p_string.dat",action="write")
   do i=1,N
    write(111,"("//int_to_str(N)//"F12.5)") p(i,:)
   end do
  close(111)
 end if
 
 !print *,"Writing p-matrix..."
 !open(111,file="p.dat",action="write")
 ! do i=1,N
 !  write(111,"(<N>F10.5)") vectors(i,1:N)
 ! end do
 !close(111)
 
 call do_clustering(p,c0,i,j,k,prom,selected_atoms,doISE,printTCL)
 !Now i contains the number of natural clusters and c0(i,:) the clusters
 print *,"The natural number of clusters is",i," cover: ",j," ratio ",prom
 !print "(<N>I4)", int(c0(i,1:N))

 !!!!

  !Prepare clusters for output
  allocate(vec(1:N))

  ! Output results
  open(111,file="clusters.hccp",action="write")
  !write name of file
  write(111,*) trim(filename)
  !write resid's 
  do i=1,N
   write(111,"(I4,$)") selected_atoms(i)%res_num
  end do
  write(111,*)
  !write chains 
  do i=1,N
   if(selected_atoms(i)%chain/=" ")then
    write(111,"(A2,$)") selected_atoms(i)%chain
   else
    write(111,"(A2,$)") "_ "
   end if
  end do
  write(111,*)
  !write SGtable
  write(111,"("//int_to_str(N_bins)//"I4)") SGtable
  !write clusters
  do i=1,N
   if(c0(i,1)>0)then !If this number of clusters was observed
    vec(1:N) = int(c0(i,1:N))
    vec =  mark_sequencially(vec,i)
    write(111,"(I5,"//int_to_str(N)//"I5)") i,vec(1:N)
   end if
  end do
  close(111)


  if(force2domains=="yes")then
   !if(i>2)then
   ! vec(1:N) = int(c0(i,1:N))
   ! vec =  mark_sequencially(vec,i)
   !else
    vec(1:N) = int(c0(2,1:N))
    vec =  mark_sequencially(vec,2)
   !end if
  else
   vec(1:N) = int(c0(i,1:N))
   vec =  mark_sequencially(vec,i)
  end if

  if(printVEC=="yes")then
   open(111,file="code_string.dat",action="write",access="append")
   !if(i>9) write(111,*) "# More than 9 clusters!"
   write(111,"(A,"//int_to_str(N)//"I1)") trim(filename)//" ",vec(1:N)
   close(111)
  end if


  if(printDB=="yes")then
   !Initialize assignment file
   open(111,file="assignment.dat",action="write",access="append")
    write(111,*) "%"
    if(chain==" ") chain = "_"
    write(111,"(A,A2,I5)") trim(filename),chain,i
   close(111)
  end if

  if(force2domains=="yes")then
   !if(i>2)then
   ! call corr_and_energy(vec,atoms,p,i,properties)
   !else
    call corr_and_energy(vec,atoms,p,2,properties)
   !end if
  else
   call corr_and_energy(vec,atoms,p,i,properties)
  end if

  print *,"     Domain correlation: ",properties%dom_corr
  print *,"Interdomain correlation: ",properties%interdom_corr
  print *,"     Interdomain energy: ",properties%interdom_en
  print *,"           Total energy: ",properties%all_en

  !print *,properties%dom_str
  if(printDB=="yes")then
   open(111,file="database.dat",action="write",access="append")
    write(111,"(A,A2,3I5,F12.6,I5,11F12.6)") trim(filename),chain,N,i,j,prom,k,&
     properties%dom_corr,properties%interdom_corr,properties%dom_en,properties%interdom_en,properties%all_en, &
     properties%dom_size1,properties%dom_size2,properties%seg_num1,properties%seg_num2,properties%dom_corr1,properties%dom_corr2
   close(111)
  end if

  deallocate(vec)
 !!!!

 deallocate(vectors,values,c0,p,beta,b1)
 goto 333
!==============================================

 !call b_factors_ANM(c0,selected_atoms,beta)

 deallocate(vectors,values,pot,c0,beta,b1,p)
333 end program
