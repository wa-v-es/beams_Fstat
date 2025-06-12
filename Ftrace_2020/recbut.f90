
! Code converted using TO_F90 by Alan Miller
! Date: 2016-06-03  Time: 15:31:34

SUBROUTINE recbut(xr,n,delr,nrdr,lfr,hfr,yr)

!     Routine to compute the recursive filter coefficients of a NRDR pole
!     Butterworth bandpass filter, with corner frequencies at LF and HF.
!     Apply filter to data in array XR (with N points and sampling rate DELR)
!     Return filtered data as array YR.

!     The filter is normalised to a gain of one in the pass band. The theory
!     is given in "Time sequence analysis in geophysics" by E. R.
!     Kanasewich (pages 181-199).

!     See AG/355 (Douglas 1993) for further details.



REAL, INTENT(IN OUT)                     :: xr(n)
INTEGER, INTENT(IN)                      :: n
REAL, INTENT(IN OUT)                     :: delr
INTEGER, INTENT(IN)                      :: nrdr
REAL, INTENT(IN OUT)                     :: lfr
REAL, INTENT(IN OUT)                     :: hfr
REAL, INTENT(OUT)                        :: yr(n)
DOUBLE PRECISION :: dela,a,b,c,h1,h2,pi,pi2,wp,alf,hf,x,an,  &
    dw,con,fo1,fo2,fo3,tbd,thet,ang2,ang
COMPLEX*16 p,pw


DIMENSION x(n),p(50),h1(50),h2(50),fo1(50),fo2(50), fo3(50)

dela = DBLE(delr)
alf = DBLE(lfr)
hf = DBLE(hfr)
DO i = 1,n
  x(i) = DBLE(xr(i))
END DO

pi=4.0D0*DATAN(1.0D0)
pi2=2.0D0*pi

tbd=2.0D0/dela
hhf=hf
hlf=alf
hf=hf*pi2
alf=alf*pi2

!     Modify the cut-off frequencies.

hf=tbd*DTAN(hf/tbd)
alf=tbd*DTAN(alf/tbd)
dw=hf-alf
wp=hf*alf
noby2=nrdr/2
ang=pi/dfloat(nrdr*2)
ang2=ang*2.0D0
con=1.0D0/dw**nrdr

!     Loop to set up poles and filter coefficients using equations
!     13.6.17 and 13.7.7.

DO  i=1,noby2
  ii=nrdr-i+1
  thet=dfloat(i-1)*ang2+ang
  pw=DCMPLX(-DCOS(thet),DSIN(thet))*dw/2.0D0
  p(i)=CDSQRT(pw**2-DCMPLX(wp,0.0D0))
  p(nrdr+i)=pw-p(i)
  p(i)=pw+p(i)
  p(ii)=DCONJG(p(i))
  b=dreal(-(p(i)+p(ii)))
  c=dreal(p(i)*p(ii))
  a=2.0D0/dela+b+c*dela/2.0D0
  h1(i)=(c*dela-4.0D0/dela)/a
  h2(i)=(2.0D0/dela-b+c*dela/2.0D0)/a
  con=con*a
  p(nrdr+ii)=DCONJG(p(nrdr+i))
  b=dreal(-(p(nrdr+i)+p(nrdr+ii)))
  c=dreal(p(nrdr+i)*p(nrdr+ii))
  a=2.0D0/dela+b+c*dela/2.0D0
  h1(noby2+i)=(c*dela-4.0D0/dela)/a
  h2(noby2+i)=(2.0D0/dela-b+c*dela/2.0D0)/a
  con=con*a
END DO
con=1.0D0/con

!     P(1) to P(NRDR) are the poles for the low frequency cut-off and
!     P(NRDR+1) to P(2*NRDR) the poles for the high frequency cut-off.
!     A,B,C are equivalent to a(sub j ),b(sub j) and c(sub j) for stage
!     J of the filter (see equation 13.7.7 of Kanasewich).
!     H1(J) and H2(J) are the filter coefficients for stage J. CON is
!     the normalising factor.


!     Apply NRDR stage recursive filter to X array. FO1(J) is output of
!     stage J at time I. FO2(J) is output at time I-1 and FO3(J) is output
!     at time I-2. The Y array is the filter output.

nlim=nrdr
DO  j=1,nlim
  fo1(j)=0.0
  fo2(j)=0.0
  fo3(j)=0.0
END DO
nlim=nrdr+1
DO  i=3,n
  fo1(1)=x(i)
  fo3(1)=x(i-2)
  DO  j=2,nlim
    fo1(j)=fo1(j-1)-fo3(j-1)-fo2(j)*h1(j-1)-fo3(j)*h2(j-1)
  END DO
  DO  j=2,nlim
    fo3(j)=fo2(j)
    fo2(j)=fo1(j)
  END DO
  yr(i)=SNGL(fo1(nlim)*con)
END DO
1     CONTINUE

RETURN
END SUBROUTINE recbut
