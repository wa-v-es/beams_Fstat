#
# Makefile for program ftrace (SAC data file version)
#
F90=gfortran
FFLAGS= -finit-local-zero -fno-automatic -O3
SACLIB=/usr/local/sac/lib/sacio.a
BINDIR=/Users/keyser/Documents/softies/SmKS-main/bin/

#Make ftrace with radpat norm, circ-wave correction, and elevation correction
CONOBJ=main_rnorm_tau_tmin_circ_elv.f90 geodesic.f90 beta_as63.f90 fstuff.f90 recbut.f90 shift.f90 dcool.f90 ctaper.f90 linreg_sub.f90
BINFILE=ftrace_rnorm_tau_tmin_circ_elv_f90

RM=/bin/rm -f

#
# Fortran programs
#

all:	${CONOBJ}
	$(F90) $(FFLAGS) -c beta_as63.f90
	$(F90) $(FFLAGS) -c geodesic.f90
	$(F90) $(FFLAGS) -o ${BINFILE} ${CONOBJ} ${SACLIB}
	cp ${BINFILE} $(BINDIR)

clean :
	$(RM) ${BINFILE} ${BINFILE_EL} ${BINFILE_CIRC} *.mod *.o
