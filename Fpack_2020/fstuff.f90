SUBROUTINE fstuff(m,spts,snr,hf,lf,npts,narr,MAX,del,DATA,semb,  &
        fst,prob)
use LogBeta
implicit none
! Code converted using TO_F90 by Alan Miller
! Date: 2016-06-03  Time: 15:31:00


INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN)                      :: spts
REAL, INTENT(IN OUT)                     :: snr
REAL, INTENT(IN OUT)                     :: hf
REAL, INTENT(IN OUT)                     :: lf
INTEGER, INTENT(IN)                      :: npts
INTEGER, INTENT(IN OUT)                  :: narr
INTEGER, INTENT(IN OUT)                  :: MAX
REAL, INTENT(IN)                         :: del
REAL, INTENT(IN)                         :: DATA(narr,MAX)
REAL, INTENT(OUT)                        :: semb(MAX)
REAL, INTENT(OUT)                        :: fst(MAX)
REAL (kind=4), INTENT(OUT)                        :: prob(MAX)
REAL :: fnn1,fprime,dum0,dum1,dum2,dum3
INTEGER :: j,jj
real (kind = 8) :: beta_log,a,b,f,tempvar,betain
!REAL :: betai

INTEGER :: smv,nwin,nn1,nn2,nc1,lambda
INTEGER :: kk

!open(3,file="probout.txt",status='replace')
! calculate semblance for a moving 2*spts+1 sample window
! and probability of the presence of a signal in band-limited white noise.
! dof nn1 = 2BT, where B = hf - lf (Hz), T = time window (sec).
!     nn2 = nn1(m-1)
! non-centrality parameter, lambda = nn1 * (snr)**2 (Blandford 1974)

nwin = 2*spts
fnn1 = 2.*(hf-lf)*nwin*del
nn1 = INT(fnn1)
IF(nn1 < 1) nn1 = 1
nn2 = nn1*(m-1)

! probability of a non-central F-distribution;
! Abramowitz and Stegun (1964) P(F'|nn1,nn2,lambda) = P(F|nc1,nn2)
! where F = (nn1 * F')/(nn1 + lambda), and
!       nc1 = (nn1 + lambda)**2/(nn1 + 2*lambda)

lambda = INT(fnn1*(snr)**2)
nc1 = (nn1 + lambda)**2/(nn1 + 2*lambda)
!WRITE(0,*) 'SNR sqrt(signal/noise power):',snr
!WRITE(0,*) 'dof:',nn1,nn2,'  lambda:',lambda
! loop round data calculating the semblance.

DO kk = 1, npts-nwin
  smv = kk + spts
  
! calculate semblance
  
! top: dum1 = sum of square of sum
! bottom: dum2 = sum of sum of squares
  dum1 = 0.0
  dum3 = 0.0
  DO jj = kk,kk+nwin
    dum0 = 0.0
    dum2 = 0.0
    DO j = 1,m
      dum0 = dum0 + DATA(j,jj)
      dum2 = dum2 + DATA(j,jj)**2
    END DO
    dum1 = dum1 + dum0**2
    dum3 = dum3 + dum2
  END DO
! semblance
  semb(smv) = dum1 / (FLOAT(m) * dum3)
! F-statistic F = S(m-1) / (1-S)
  fst(smv) = semb(smv) * FLOAT(m-1) / (1. - semb(smv))
! probability of a non-central F-distribution.
  
  fprime = (FLOAT(nn1)*fst(smv))/FLOAT(nn1 + lambda)


  a=0.5*nn2
  b=0.5*nc1
  f=(nn2/(nn2+nc1*fprime))

  beta_log= betaln(a,b)

!Convert Real*8 back into Real*4 for writing out ('wsac' can only write Real*4)  - MODIFIED BY DAN FROST (26.10.16)
  prob(smv) = real((1 - betain(f,a,b,beta_log)),4)


END DO

RETURN
END SUBROUTINE fstuff
