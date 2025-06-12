
! Code converted using TO_F90 by Alan Miller
! Date: 2016-06-03  Time: 15:31:38

! Subroutine SHIFT(Npts,Array,Npts2,Time,Buffer)

!***************************************************************************

! This subroutine uses a Fourier transform to time shift an array of
! data. The input data is transformed in and out of the frequency
! domain using the Dcool subroutine and the output data replaces the
! input data in Array. This subroutine is based on subroutine SFCH
! by John Young, Blacknest.

!***************************************************************************

! Programmer: Howard Rainford
! Date:       14 June 1995

!***************************************************************************

!       Called by:
!  Time_Shift

!       Usage:  Call Shift(Npts,Array,Npts2,Time,Buffer)

!        Arguments:
!  Npts : i*4    Number of points in Array
!  Array : real*4 Input data array (returned with time shift)
!  Npts2   : i*4    Number corresponding to the next highest
!     power of two above npts
!  Time : real*4 Time in seconds by which trace is to be
!     shifted
!  Buffer :cmplx*16 Complex work array

! Subroutine and Function Calls:
!  Dcool  Performs Fourier and Inverse Fourier
!    transforms on a complex data array

!***************************************************************************

! Variables Used:
! dnorm           :  Normalising constant (complex Npts)
! i,ii  : Loop counters
! nhalf           :  Npts/2 + 1
! n2              :  Integer power of two corresponding to Npts2
! omega           :  Complex calculation constant
! t               :  Complex time shift variable
! twopi           :  2*pi

!***************************************************************************
! Modification History

! DATE       VERSION AUTHOR
! 14-JUN-1995     1.0 Howard Rainford
!  First Release

!***************************************************************************

SUBROUTINE shift(npts,array,npts2,time,buffer)

IMPLICIT NONE
!--- All variables must be explicitly declared

INTEGER*4, INTENT(IN)                    :: npts
REAL*4, INTENT(OUT)                      :: array(npts)
INTEGER*4, INTENT(IN)                    :: npts2
REAL*4, INTENT(IN OUT)                   :: time
COMPLEX*16, INTENT(OUT)                  :: buffer(npts2)


!--- Arguments






!--- Local variables
INTEGER*4  i, ii, n2, nhalf
REAL*8     twopi
COMPLEX*16 omega, t, dnorm


!------- Routine Start -----------------------------------------------


!--- Get some constants  &
nhalf = npts2/2 + 1
dnorm = DCMPLX(npts2)
twopi = 8.0D0 * DATAN(1.0D0)
n2    = nint(ALOG10(FLOAT(npts2))/ALOG10(2.0))


!--- Fill the real part of buffer with the input data and then pad to
!--- npts2 with zeroes
DO i = 1, npts
  buffer(i) = DCMPLX(array(i))
END DO

DO i = npts+1, npts2
  buffer(i) = DCMPLX(0.0)
END DO


!--- Perform a Fourier transform on buffer
call dcool(n2,buffer,-1.0)


!--- Get some more constants
omega = DCMPLX(0.0D0,DBLE(time)*twopi) / dnorm
omega = CDEXP(omega)
t     = DCMPLX(1.0,0.0)


!--- Do the time shift in the frequency domain
buffer(1) = buffer(1) / dnorm
DO i = 2, nhalf
  ii = npts2 + 2 - i
  t  = t * omega
  buffer(i)  = buffer(i) * t / dnorm
  buffer(ii) = DCONJG(buffer(i))
END DO


!--- Perform an inverse Fourier transform on buffer
!call myFFT(buffer,n2,-1.0)

!call dcool(n2,buffer,-1.0)
CALL dcool(n2,buffer,1.0)


!--- Output array is the real part of buffer
DO i = 1, npts
  array(i) = dreal(buffer(i))
END DO

RETURN

END SUBROUTINE shift
