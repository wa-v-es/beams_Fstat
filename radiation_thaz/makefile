#Makefile for program radiathion_thaz

#Variables
#---------
F90=gfortran
FFLAGS=-O3 
RM=/bin/rm -f
BINDIR=/Users/keyser/Documents/softies/SmKS-main/bin/

all : main

#Compile Source-code and link.
radiation_thaz : radiation_thaz.f
	$(F90) $(FFLAGS) radiation_thaz.f -o radiation_thaz

#Copy executable and man page to appropriate directories.
main : radiation_thaz
	cp radiation_thaz $(BINDIR)

clean :
	  $(RM) radiation_thaz

