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
REAL :: DATA,x,y,tau,array,semb,fst,beg,del !,prob
real (kind=4) :: prob
REAL :: lat,lon,elv
DOUBLE PRECISION :: dist,az,dtokm
real (kind=4) :: clat, clon,elat,elon
DOUBLE PRECISION :: dclat,dclon,dlat,dlon,ddist,daz,da1,dd1
DOUBLE PRECISION :: cdist,az_cgeod,cdaz,dist_rel
DOUBLE PRECISION :: a,b,chi2
REAL :: baz,pvel,azb,skx,sky
REAL :: dtor,twopi
REAL :: an2,snr,sta,lta,det,fnn1,fprime
REAL :: ang,dsq
REAL :: siga,sigb,xx,q
REAL :: xir,filt,lf,hf
REAL :: time,dummy
REAL :: tmin,elv_vel
CHARACTER :: DO_ELV, DO_CIRC, DO_RADPAT
real (kind=8) :: a_geod,f_geod
real (kind=8) :: az_geod,dum1,dum2,dum3,dum4,dum5,dum6,dum7
integer :: dumo
integer :: masko
real maxnum,maxnum2,minnum,minnum2
integer radpat_updo
!integer, parameter :: sp = kind(1e0)

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
CHARACTER (LEN=6) :: kstnm_var(narr)
CHARACTER (LEN=80) :: kname(narr)
DIMENSION DATA(narr,MAX)
DIMENSION x(narr),y(narr),lat(narr),lon(narr),elv(narr),dist(narr),az(narr),az_geod(narr)
DIMENSION dlat(narr),dlon(narr),ddist(narr),daz(narr)
DIMENSION tau(narr),array(MAX),semb(MAX),xir(MAX),filt(MAX)
DIMENSION time(MAX),dummy(MAX),prob(MAX),fst(MAX),det(MAX)
DIMENSION w(2*MAX)
DOUBLE PRECISION :: elat_out,elon_out
LOGICAL :: arcmod

! linux g77 features
ifc=system('/bin/rm -f beam.sac filt.sac semb.sac xf.sac xp.sac')
kcmpnm=' '
! constants
twopi = 8.00 * ATAN(1.0)
dtor = twopi / 360.0

! input program control information and data:
! read in baz, pvel and filem

open(4,file='tau.output')
READ(*,*) m,baz,pvel
READ(*,*) clat,clon
READ(*,*) elat,elon
DO i = 1,m
  READ(*,99) kname(i)
END DO
99   FORMAT(a80)
READ(*,*) spts,lf,hf,npol,snr
read(*,*) tmin
read (*,*) DO_CIRC, DO_ELV, elv_vel, DO_RADPAT, radpat_updo

! Central dist
call invers(a_geod,f_geod,dble(elat),dble(elon),dble(clat),dble(clon),cdist,az_cgeod,cdaz&
     ,masko,dum1,dum2,dum3,dum4,dum5)
write(*,*) "cdist",cdist

!New elat
arcmod=.false.
call direct(a_geod,f_geod,dble(clat),dble(clon),dble(baz),cdist*1000,arcmod&
     ,elat_out,elon_out,dum1,dumo,dum3,dum4,dum5,dum6,dum7)
elat=elat_out
elon=elon_out
! Central dist


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


IF (pvel > 0.0) THEN
  azb = baz * dtor
  skx = SIN(azb) / pvel
  sky = COS(azb) / pvel
END IF
dclat = DBLE(clat * dtor)
dclon = DBLE(clon * dtor)
!DO i = 1,m
DO i = 1,m
  CALL rsac1(kname(i),array,npts,beg,del,MAX,nerr)
  CALL getfhv('stla',lat(i),nerr)
  CALL getfhv('stlo',lon(i),nerr)
  IF ( DO_ELV == "Y" ) THEN
     CALL getfhv('stel',elv(i),nerr)
  END IF
  CALL getkhv('kstnm',kstnm,nerr)
  kstnm_var(i)=trim(kstnm)
  IF(i == 1) THEN
    CALL getnhv('nzyear',nzyear,nerr)
    CALL getnhv('nzjday',nzjday,nerr)
    CALL getnhv('nzhour',nzhour,nerr)
    CALL getnhv('nzmin',nzmin,nerr)
    CALL getnhv('nzsec',nzsec,nerr)
    CALL getnhv('nzmsec',nzmsec,nerr)
  END IF
!DEBUG  WRITE(0,*) lat(i),lon(i)
  dlat(i) = DBLE(lat(i) * dtor)
  dlon(i) = DBLE(lon(i) * dtor)
  daz(i)=0; ddist(i)=0


  !CIRCULAR WAVE FRONT
  IF ( DO_CIRC == "Y" ) THEN
  call invers(a_geod,f_geod,dble(lat(i)),dble(lon(i)),dble(elat),dble(elon),ddist(i),az_geod(i),daz(i)&
  ,masko,dum1,dum2,dum3,dum4,dum5)
  dist_rel=ddist(i)-cdist

  IF ( DO_ELV == "Y" ) THEN
     elv(i)=elv(i)/1000
     tau(i) = (dist_rel)/pvel - elv(i)/elv_vel
  ELSE
     tau(i) = (dist_rel)/pvel
  END IF
  WRITE(*,91) kstnm_var(i),tau(i)


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
  WRITE(*,92) kstnm_var(i),tau(i),x(i),y(i)

  END IF
  92 FORMAT('Station, Tau, X, Y: ',4(A,3f11.3))
  91 FORMAT('Station, Tau: ',2(A,f11.5))
  !Write stanm and time shift out to tau.output
  WRITE(4,91) kstnm_var(i),tau(i)


! remove linear trend, 1% cosine taper, time shift and load into data array
! Change beginning time to that read from input file - solves problem with beams coming out with start=0
  beg=tmin

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


! Calculate linear beam (unfiltered and un-normalised) - delay and sum data array
DO ii = 1,npts
  array(ii) = 0.0
END DO
DO i = 1,m
  DO ii = 1,npts
    array(ii) = array(ii) + DATA(i,ii)
  END DO
END DO
DO ii = 1,npts
  array(ii) = array(ii) / float(m)
END DO

!Write out linear beam (unfiltered and un-normalised)
CALL setlhv('leven',.true.,nerr)
CALL setfhv('b',beg,nerr)
CALL setfhv('delta',del,nerr)
CALL setnhv('npts',npts,nerr)
CALL setnhv('nzyear',nzyear,nerr)
CALL setnhv('nzjday',nzjday,nerr)
CALL setnhv('nzhour',nzhour,nerr)
CALL setnhv('nzmin',nzmin,nerr)
CALL setnhv('nzsec',nzsec,nerr)
CALL setnhv('nzmsec',nzmsec,nerr)
CALL setkhv('kstnm',kstnm,nerr)
CALL setfhv('stla',clat,nerr)
CALL setfhv('stlo',clon,nerr)
CALL setkhv('kcmpnm',kcmpnm,nerr)
kuser1='BEAM'
CALL setkhv('kuser1',kuser1,nerr)
CALL wsac0('beam.sac',dummy,array,nerr)


!Obsolete ! calculate filtered beam (Simply filters existing beam of unfiltered data, "array")
!Obsolete ! Filter stacked data: corner frequencies at lf and hf with npol butterworth filter.
!Obsolete DO j = 1,npts
!Obsolete   jj = npts - j + 1
!Obsolete   xir(j) = array(jj)
!Obsolete END DO
!Obsolete CALL recbut(xir,npts,del,npol,lf,hf,filt)
!Obsolete DO j = 1,npts
!Obsolete   jj = npts - j + 1
!Obsolete   xir(j) = filt(jj)
!Obsolete END DO
!Obsolete CALL recbut(xir,npts,del,npol,lf,hf,filt)
!Obsolete kuser1='FILTERED'
!Obsolete call setkhv('kuser1',kuser1,nerr)
!Obsolete call wsac0('filt.sac',dummy,filt,nerr)



! Calculate filtered data, and then stack to beam
! Filter all data: corner frequencies at lf and hf with npol butterworth filter.
!
!Loop over traces
do i = 1,m
! two-pass filtering (load trace in reverse order)
   do j = 1,npts
      jj = npts - j + 1
      xir(j) = data(i,jj)
   enddo
   call recbut(xir,npts,del,npol,lf,hf,filt)
   do j = 1,npts
      jj = npts - j + 1
      xir(j) = filt(jj)
   enddo
   call recbut(xir,npts,del,npol,lf,hf,filt)

! normalise filtered trace
   maxnum = abs(maxval(filt))
   minnum = abs(minval(filt))
   !Normalise to amplitude of rad pat (single rad pat for array centre)
   IF ( DO_RADPAT == "Y" ) THEN
      write(*,*) "Radpat Y",radpat_updo
      if (radpat_updo .eq. 1) then
         filt = filt/maxnum
      else if  (radpat_updo .eq. -1) then
         filt = filt/(abs(minnum))
      endif
      !Normalise to greatest absolute amplitude
   ELSE IF ( DO_RADPAT == "M" ) THEN
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
   maxnum2 = maxval(filt)
   minnum2 = minval(filt)

   print *,"m=",i,"max",maxnum,"maxnorm=",maxnum2,"minnorm=",minnum2
! load into data array
   do j = 1,npts
      data(i,j) = filt(j)

   enddo
enddo



!-----------------------------------------------------
! calculate beam of normalised, filtered data (D Frost 18.9.2015)
!
      do ii = 1,npts
        array(ii) = 0.0
      enddo
      do i = 1,m
         do ii = 1,npts
            array(ii) = array(ii) + data(i,ii)
         enddo
      enddo
      do ii = 1,npts
        array(ii) = array(ii) / float(m)
      enddo
          kuser1='BeamOfFiltData'
          call setkhv('kuser1',kuser1,nerr)
!      call wsac0('beamfilt.sac',dummy,array,nerr)
      call wsac0('filt.sac',dummy,array,nerr)
!-----------------------------------------------------

CALL fstuff(m,spts,snr,hf,lf,npts,narr,MAX,del,DATA,semb,fst,prob)

! cosine taper and write out semblance, F, probability


nn = 2*spts+1
CALL ctaper(npts,semb,nn)
kuser1='SEMBLANCE'
CALL setkhv('kuser1',kuser1,nerr)
CALL wsac0('semb.sac',dummy,semb,nerr)
kuser1='F TRACE'
CALL setkhv('kuser1',kuser1,nerr)
CALL wsac0('xf.sac',dummy,fst,nerr)
kuser1='P(F)'
CALL setkhv('kuser1',kuser1,nerr)
CALL wsac0('xp.sac',dummy,prob,nerr)

STOP
END PROGRAM ftrace
