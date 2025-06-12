
! Code converted using TO_F90 by Alan Miller
! Date: 2016-06-03  Time: 15:31:22

SUBROUTINE ctaper(asize,aa,nn)


INTEGER, INTENT(IN OUT)                  :: asize
REAL, INTENT(OUT)                        :: aa(asize)
INTEGER, INTENT(IN OUT)                  :: nn
INTEGER :: i
REAL :: pi, ang

pi = 4.00 * ATAN(1.0)
ang=pi/(2.*nn)
DO  i=1,nn
  aa(nn-i+1)=aa(nn-i+1)*((COS((1-i)*ang))**2)
  aa(asize-nn+i)=aa(asize-nn+i)*((COS((1-i)*ang))**2)
END DO

RETURN
END SUBROUTINE ctaper
