PROGRAM ftrace

! Code converted using TO_F90 by Alan Miller
! Date: 2016-06-03  Time: 14:49:12

!   David Bowers 29/6/2005
!   Modified by SP/DB 4/1/06 to make number of poles of filter always
!   even - by adding one to it if it is odd, and warning user.
! Program to calculate and apply a time shift given a backazimuth (BAZ),
! phase velocity (PVEL) of a series of M array waveforms (up to narr=55).
! x and y offsets are calculated from station latitudes and longitudes
! relative to the array centre (CLAT, CLON).

! Bandpass filter data (using a recursive Butterworth filter)
! with corner frequencies at LF and HF, and NPOL order (must be even).
! (two-pass acausal filtering).

! Calculates F-trace using a moving window 2*SPTS+1 long.

! Probability (using the F-statistic) that a signal is present above
! bandlimited noise, with SNR equal to sqrt(signal/noise power)
! (see Douze and Laster, 1979. Geophysics. 44, pp.1999 and
! Blandford, 1974. Geophysics. 39, pp.633)

! f77 -fast -o ftrace ftrace.f /usr/remote/apps/sac/lib/sac.a
! /usr/remote/apps/lib/librecipes_f.a

! cobbled together by Dave Bowers

! input: SAC format time series files with STLA and STLO defined.

! control information is read from standard input
! (NB. input filenames are character*40)

! --- input file ---
! m baz pvel
! clat clon
! file1
! file2
!   .    .  .
! filem                                      (format a60)
! spts lf hf npol snr
! ---


! output: delay-and-sum 'beam.sac'
!         filtered delay-and-sum 'filt.sac'
!         F-detector as a function of time 'xf.sac'
!         Probability of a signal with SNR 'xp.sac'

! Uses:  shift
!        dcool
!        rsac1  etc. (apps/sac/lib/sac.a)
!        fit         (apps/lib/librecipes_f.a)
!        ctaper
!        recbut
!        george
!        fstuff
!        betai       (/apps/lib/librecipes_f.a)

use geodesic

IMPLICIT NONE
REAL :: DATA,x,y,tau,array,semb,fst,beg,del,array_all !,prob
real (kind=4) :: prob
REAL :: lat,lon,elv
DOUBLE PRECISION :: dist,az,dtokm
real (kind=4) :: clat, clon,elat,elon
DOUBLE PRECISION :: dclat,dclon,dlat,dlon,ddist,daz,da1,dd1
DOUBLE PRECISION :: cdist,az_cgeod,cdaz,dist_rel
DOUBLE PRECISION :: a,b,chi2
REAL :: pvel,azb,skx,sky
integer rel_baz,rel_baz_min,rel_baz_inc,rel_baz_max,gcp_RD,baz
integer uloop_val,uloop_min,uloop_max,uloop_inc
real u_min,u_inc,u_max,u_val
REAL :: dtor,twopi
REAL :: an2,snr,sta,lta,det,fnn1,fprime
REAL :: ang,dsq
REAL :: siga,sigb,xx,q
REAL :: xir,filt,lf,hf
REAL :: time,dummy,tmin,elv_vel
CHARACTER :: DO_ELV, DO_CIRC, DO_RADPAT
real (kind=8) :: a_geod,f_geod
real (kind=8) :: az_geod,dum1,dum2,dum3,dum4,dum5,dum6,dum7
integer :: dumo
integer :: masko,radpat_updo
real maxnum,minnum
parameter(a_geod=6356774.719,f_geod=1/298.25) !Australia ellipsoid used in george! 
!parameter(a_geod=6356752.3142,f_geod=1/298.257223563)
parameter(masko=1)
COMPLEX*16 w
INTEGER :: nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec
INTEGER :: npts,nfil,npol,nhalf
INTEGER :: i,ii,j,jj,k,kk,l,ll,m,mm,n,nn,n2,npts2,nerr,MAX,narr
INTEGER :: spts,smv,nwin,nn1,nn2,nc1,lambda
INTEGER :: llta,rlta,df
INTEGER :: ifc, system
PARAMETER(narr=150,MAX=131072,nfil=2048)
PARAMETER(dtokm=111.19504)
CHARACTER (LEN=8) :: kstnm,kcmpnm
CHARACTER (LEN=8) :: kuser1
CHARACTER (LEN=80) :: kname(narr),fileout,form
DIMENSION DATA(narr,MAX)
DIMENSION x(narr),y(narr),lat(narr),lon(narr),elv(narr),dist(narr),az(narr),az_geod(narr)
DIMENSION dlat(narr),dlon(narr),ddist(narr),daz(narr)
DIMENSION tau(narr),array(MAX),semb(MAX),xir(MAX),filt(MAX),array_all(narr,MAX)
DIMENSION time(MAX),dummy(MAX),prob(MAX),fst(MAX),det(MAX)
DIMENSION w(2*MAX)
DOUBLE PRECISION :: elat_out,elon_out
LOGICAL :: arcmod

! linux g77 features
kcmpnm=' '
! constants
twopi = 8.00 * ATAN(1.0)
dtor = twopi / 360.0

! input program control information and data:
! read in rel_baz, pvel and filem

read(*,*) m
read(*,*) rel_baz_min,rel_baz_inc,rel_baz_max
read(*,*) u_min,u_inc,u_max
write(*,*) 'rel_baz_min',rel_baz_min,'rel_baz_inc',rel_baz_inc,'rel_baz_max',rel_baz_max
write(*,*) 'u_min',u_min,'u_inc',u_inc,'u_max',u_max
READ(*,*) clat,clon
READ(*,*) elat,elon
DO i = 1,m
  READ(*,99) kname(i)
END DO
99   FORMAT(a80)
READ(*,*) spts,lf,hf,npol,snr
read(*,99) fileout
read(*,*) tmin,gcp_RD
read (*,*) DO_CIRC, DO_ELV, elv_vel, DO_RADPAT, radpat_updo


! Central dist
call invers(a_geod,f_geod,dble(elat),dble(elon),dble(clat),dble(clon),cdist,az_cgeod,cdaz&
     ,masko,dum1,dum2,dum3,dum4,dum5)
!!  write(*,*) "invers in",elat,elon,clat,clon,"out",cdist,az_cgeod
!!write(*,*) "cdist",cdist
! Central dist


!     Prepare output file
!     relbaz, u_val, time, filt, fst baz
      form="(I5 F6.2 F8.2 F20.10 F20.10 I5)"

      open(3,file=fileout,status='replace')


!     check that number of poles of filter is greater than zero; if
!     not, set it to default value of two; else check that it is even:
!     if not, add one to it. Warn user in either case - SP/DB 4/1/2006
IF (npol <= 0) THEN
  WRITE (6, 1030) npol
  1030   FORMAT ("WARNING: Filter number of poles, ", i4,  &
      ", is zero or less:")
  npol = 2
  WRITE (6, 1040) npol
  1040   FORMAT ("WARNING: it has been set to ", i4, ".")
ELSE IF (MOD(npol,2) /= 0) THEN
  WRITE (6, 1050) npol
  1050   FORMAT ("WARNING: Filter number of poles, ",i4,", is not even:")
  npol = npol + 1
  WRITE (6, 1060) npol
  1060   FORMAT ("WARNING: it has been increased to ", i4, ".")
END IF
! Read station lat and long, calculate offset (relative to CLAT, CLON).
! Calculate time delays tau(i) due to array element offset x(i),y(i)

!----
!READ IN DATA
DO i = 1,m
   CALL rsac1(kname(i),array,npts,beg,del,MAX,nerr)
   DO j = 1,npts
      array_all(i,j)=array(j)
   END DO

  CALL getfhv('stla',lat(i),nerr)
  CALL getfhv('stlo',lon(i),nerr)
  IF ( DO_ELV == "Y" ) THEN
     CALL getfhv('stel',elv(i),nerr)
  END IF
  IF(i == 1) THEN
    CALL getnhv('nzyear',nzyear,nerr)
    CALL getnhv('nzjday',nzjday,nerr)
    CALL getnhv('nzhour',nzhour,nerr)
    CALL getnhv('nzmin',nzmin,nerr)
    CALL getnhv('nzsec',nzsec,nerr)
    CALL getnhv('nzmsec',nzmsec,nerr)
    CALL getkhv('kstnm',kstnm,nerr)
  END IF
!DEBUG  WRITE(0,*) lat(i),lon(i)
  dlat(i) = DBLE(lat(i) * dtor)
  dlon(i) = DBLE(lon(i) * dtor)
  daz(i)=0; ddist(i)=0
END DO
!----


uloop_max=((u_max-u_min)/u_inc)
uloop_min=0
uloop_inc=1

elat=0
elon=0
do 101 rel_baz = rel_baz_min,rel_baz_max,rel_baz_inc

!            if ((baz-gcp_RD) .le. -180) then
!               relbaz=((baz-gcp_RD)+360)
!            else if ((baz-gcp_RD) .gt. 180) then
!               relbaz=((baz-gcp_RD)-360)
!            else
!               relbaz=(baz-gcp_RD)
!            end if

            if ((gcp_RD+rel_baz) .le. -180) then
               baz=(360+(gcp_RD+rel_baz))
            else if ((gcp_RD+rel_baz) .gt. 180) then
               baz=((gcp_RD+rel_baz)-360)
            else
               baz=(gcp_RD+rel_baz)
            end if
!               write(*,*) 'baz',baz,"cdist",cdist
! Central dist
arcmod=.false.
call direct(a_geod,f_geod,dble(clat),dble(clon),dble(baz),cdist*1000,arcmod&
     ,elat_out,elon_out,dum1,dumo,dum3,dum4,dum5,dum6,dum7)
!call direct(a_geod,f_geod,dble(clat),dble(clon),dble(baz),cdist,arcmod&
!     ,dble(elat),dble(elon),dum1,dumo,dum3,dum4,dum5,dum6,dum7)
elat=elat_out
elon=elon_out
!!write(*,*) "direct in",clat,clon,baz,cdist,"out",elat,elon
! Central dist
!!write(*,*) "elat",elat,"elon",elon
!!write(*,*) "elat_out",elat_out,"elon_out",elon_out
!      do 101 u_val = u_min,u_max,u_inc

do 102 uloop_val = uloop_min,uloop_max,uloop_inc
   u_val=u_min+(u_inc*uloop_val)
   if (u_val.eq.0.0) then
      pvel=111.1/0.001
   else 
      pvel=111.1/u_val
   endif

!   write(*,*) 'pvel',pvel,'u_val',u_val
!   write(*,*) 'rel_baz_min',rel_baz_min,'rel_baz_inc',rel_baz_inc,'rel_baz_max',rel_baz_max


IF (pvel > 0.0) THEN
  azb = baz * dtor
  skx = SIN(azb) / pvel
  sky = COS(azb) / pvel
END IF
dclat = DBLE(clat * dtor)
dclon = DBLE(clon * dtor)
!DO i = 1,m
DO i = 1,m

   DO j = 1,npts
      array(j)=array_all(i,j)
   END DO
   
  
  !CIRCULAR WAVE FRONT
  IF ( DO_CIRC == "Y" ) THEN
  call invers(a_geod,f_geod,dble(lat(i)),dble(lon(i)),dble(elat),dble(elon),ddist(i),az_geod(i),daz(i)&
       ,masko,dum1,dum2,dum3,dum4,dum5)
!!  write(*,*) "invers in",lat(i),lon(i),elat,elon,"out",ddist(i),az_geod(i)
  dist_rel=ddist(i)-cdist
!!  write(*,*) ddist(i),cdist,dist_rel
  IF ( DO_ELV == "Y" ) THEN 
     elv(i)=elv(i)/1000
     tau(i) = (dist_rel)/pvel - elv(i)/elv_vel
  ELSE
     tau(i) = (dist_rel)/pvel
  END IF
!  write(*,*) "dist_rel",dist_rel
!  WRITE(*,*) kstnm,tau(i),dist_rel
!  91   FORMAT('Station, Tau,: ',3(A,2f11.3))


  !PLANE WAVE FRONT
  ELSE IF ( DO_CIRC == "N" ) THEN
  call invers(a_geod,f_geod,dble(clat),dble(clon),dble(lat(i)),dble(lon(i)),ddist(i),az_geod(i),daz(i)&
       ,masko,dum1,dum2,dum3,dum4,dum5)
  ang = SNGL(daz(i))
  dist(i) = SNGL(ddist(i))
  x(i) = dist(i) * SIN(ang)
  y(i) = dist(i) * COS(ang)

  IF ( DO_ELV == "Y" ) THEN 
     elv(i)=elv(i)/1000
     tau(i) = ((-x(i) * skx) - (y(i) * sky)) - elv(i)/elv_vel
  ELSE
     tau(i) = ((-x(i) * skx) - (y(i) * sky))
  END IF
!  WRITE(*,92) kstnm,tau(i),x(i),y(i)

  END IF
!  92 FORMAT('Station, Tau, X, Y: ',4(A,3f11.3))
!  91 FORMAT('Station, Tau: ',2(A,1f11.5))


  

! remove linear trend, 1% cosine taper, time shift and load into data array
  DO ii = 1,npts
    time(ii)= beg + del * real(ii-1)
  END DO
  CALL linreg_sub(time,array,npts,a,b,chi2)
!DBEUG  write(*,*) 'b',b,'a',a,'chi2',chi2

  xx = beg
  DO j = 1,npts
    array(j) = array(j) - a - b * xx
!           write(*,*) "i",i,"j",j,"ar",array(j),"b",b,"a",a,"xx",xx
    xx = xx + del
  END DO
  nn = nint(real(npts)*0.05)


  CALL ctaper(npts,array,nn)
! dimensionless units for tau
  tau(i) = tau(i) / del

  
! calculation of npts2 for shift routine
  
  an2 = ALOG10(real(npts-1))/ALOG10(2.0)
  n2 = INT(an2) + 1
  npts2 = 2**n2
  nhalf = npts2/2 + 1
 
!!DEBUG - Write out array 
!DO ii = 1,npts
!write(*,*) 'array',array(ii)
!END DO

  CALL shift(npts,array,npts2,tau(i),w)

!!DEBUG - Write out shifted array 
!DO ii = 1,npts
!write(*,*) 'array',array(ii)
!END DO


! load time shifted channels into data array
  DO ii = 1,npts
    DATA(i,ii) = array(ii)
  END DO
END DO
!91   FORMAT('Tau, X, Y: ',3(2X,f11.5))

! delay and sum data array

!DO ii = 1,npts
!  array(ii) = 0.0
!END DO
!DO i = 1,m
!  DO ii = 1,npts
!    array(ii) = array(ii) + DATA(i,ii)
!  END DO
!END DO
!DO ii = 1,npts
!  array(ii) = array(ii) / float(m)
!END DO


! Calculate filtered data, and then stack to beam
! Filter all data: corner frequencies at lf and hf with npol butterworth filter.
DO i = 1,m
! two-pass filtering (load data in reverse order)
  DO j = 1,npts
    jj = npts - j + 1
    xir(j) = DATA(i,jj)
  END DO
  CALL recbut(xir,npts,del,npol,lf,hf,filt)
  DO j = 1,npts
    jj = npts - j + 1
    xir(j) = filt(jj)
  END DO
  CALL recbut(xir,npts,del,npol,lf,hf,filt)

! normalise filtered trace
   maxnum = abs(maxval(filt))
   minnum = abs(minval(filt))
   !Normalise to amplitude of rad pat (single rad pat for array centre)
   IF ( DO_RADPAT == "Y" ) THEN
      if (radpat_updo .eq. 1) then
         filt = filt/maxnum
      else if  (radpat_updo .eq. -1) then
         filt = filt/(abs(minnum))
      endif
      !Normalise to greatest absolute amplitude
   ELSE  IF ( DO_RADPAT == "M" ) THEN
      if (maxnum .gt. minnum) then
         maxnum=maxnum
      else
         maxnum=minnum
      endif
      filt = filt/maxnum
      !No normalisation
   ELSE IF ( DO_RADPAT == "N" ) THEN
      filt=filt
   ENDIF


! load into data array
  DO j = 1,npts
    DATA(i,j) = filt(j)
  END DO
END DO

! calculate filtered beam
DO ii = 1,npts
  filt(ii) = 0.0
END DO
DO i = 1,m
  DO ii = 1,npts
    filt(ii) = filt(ii) + DATA(i,ii)
  END DO
END DO
DO ii = 1,npts
  filt(ii) = filt(ii) / float(m)
END DO

!DO j = 1,npts
!  jj = npts - j + 1
!  xir(j) = array(jj)
!END DO
!CALL recbut(xir,npts,del,npol,lf,hf,filt)
!DO j = 1,npts
!  jj = npts - j + 1
!  xir(j) = filt(jj)
!END DO
!CALL recbut(xir,npts,del,npol,lf,hf,filt)
!kuser1='FILTERED'
!CALL setkhv('kuser1',kuser1,nerr)
!CALL wsac0('filt.sac',dummy,filt,nerr)


CALL fstuff(m,spts,snr,hf,lf,npts,narr,MAX,del,DATA,semb,fst,prob)

! cosine taper and write out semblance, F, probability

write(3,form) (rel_baz, u_val, ((j-1)*del)+tmin, fst(j), filt(j), baz, j=1, npts, 10)  

 102  continue

 101  continue

STOP
END PROGRAM ftrace
