#Makefile for program sac2xy

#Variables
#---------
F90=gfortran
FFLAGS=-O3 
RM=/bin/rm -f
BINDIR=/Users/keyser/Documents/softies/SmKS-main/bin/

all : main

#Compile modules.
mod_sac_io.o : mod_sac_io.f90
	$(F90) $(FFLAGS) -c mod_sac_io.f90

#Compile Source-code and link.
sachead : sachead.f90
	$(F90) $(FFLAGS) sachead.f90 -o sachead ./mod_sac_io.o

#Copy executable and man page to appropriate directories.
main : mod_sac_io.o sachead
	cp sachead $(BINDIR)

clean :
	  $(RM) sachead mod_sac_io.o sac_i_o.mod 

