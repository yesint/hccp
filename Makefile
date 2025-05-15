
FC = gfortran

COPT = -O3 -c  
OOPT = -o

hccp: comm.o nrtype.o nr.o nrutil.o eigen.o pdb.o functions.o hccp.o input_parser_util.o input_parser2.o hccp_classic.o 
	$(FC) $(OOPT) hccp.exe comm.o nrtype.o nr.o nrutil.o eigen.o pdb.o functions.o hccp.o input_parser_util.o input_parser2.o hccp_classic.o
hccp_classic.o:  hccp_classic.f90
	$(FC) $(COPT) hccp_classic.f90
input_parser2.o:  input_parser2.f90
	$(FC) $(COPT) input_parser2.f90
input_parser_util.o:  input_parser_util.f90
	$(FC) $(COPT) input_parser_util.f90
hccp.o:  hccp.f90
	$(FC) $(COPT) hccp.f90
functions.o:  functions.f90
	$(FC) $(COPT) functions.f90
pdb.o:  pdb.f90
	$(FC) $(COPT)  pdb.f90
eigen.o:  eigen.f90
	$(FC) $(COPT)  eigen.f90
nrutil.o: nrutil.f90  
	$(FC) $(COPT)  nrutil.f90 
nrtype.o:  nrtype.f90
	$(FC) $(COPT)  nrtype.f90 
nr.o:   nr.f90
	$(FC) $(COPT)  nr.f90 
comm.o: comm.f90
	$(FC) $(COPT)  comm.f90 
clean:
	rm *.o *.mod *.exe

