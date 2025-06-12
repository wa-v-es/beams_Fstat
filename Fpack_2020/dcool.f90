
! Code converted using TO_F90 by Alan Miller
! Date: 2016-06-06  Time: 22:48:58

SUBROUTINE dcool(n,xx,signi)

!     THIS SUBROUTINE WAS PROGRAMMED BY I.MACLEOD, DEPT. OF
!     ENGINEERING PHYSICS,A.N.U. AND HAS BORROWED FROM D. MCCOWAN'S
!     COOL AND IBM'S HARM.

!     DOUBLE PRECISION VERSION MODIFIED BY J.B.YOUNG FOR THE 360/75.



INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(OUT)            :: xx(2)
REAL, INTENT(IN OUT)                     :: signi
DIMENSION w(14), nbit(20),jnt(20)



DOUBLE PRECISION :: root2,pi2,con1,arg,w,cssqa,cssq2a,cssq3a
DOUBLE PRECISION :: a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i,a4r,a4i,a5r,a5i,a6r,a6i,a7r,a7i 
DOUBLE PRECISION :: a8r,a8i,xk0wr,xk0wi,xk1wr,xk1wi,xk2wr,xk2wi,xk3wr,xk3wi
DOUBLE PRECISION :: xk4wr,xk4wi,xk5wr,xk5wi,xk6wr,xk6wi,xk7wr,xk7wi,holdr,holdi


INTEGER :: offset

!nx=0
DATA nx/0/


IF(nx > 0)GO TO 100
root2=DSQRT(2.0D0)
pi2=8.0D0*DATAN(1.0D0)

100  nx=2**n
nx2=nx+nx
nx2ls1=nx2-1
nx2ls2=nx2-2
nxon8=nx/8
nxon4=nxon8+nxon8
nxon2=nxon4+nxon4
con1=pi2/dfloat(nx)

IF(signi > 0.0)GO TO 120

DO 110 k=1,nx2ls1,2
  xx(k+1)=-xx(k+1)
110 CONTINUE

120 DO 130 k=1,n
  jnt(k)=2**(n-k)
130 CONTINUE

lstart=n-n/3*3+1
IF(lstart == 1)GO TO 200
IF(lstart == 2)GO TO 150
lblok2=nxon2
l2blok=lblok2-1

DO 140 k0=1,l2blok,2
  k1=k0+lblok2
  k2=k1+lblok2
  k3=k2+lblok2
  a0r=xx(k0)+xx(k2)
  a0i=xx(k0+1)+xx(k2+1)
  a1r=xx(k0)-xx(k2)
  a1i=xx(k0+1)-xx(k2+1)
  a2r=xx(k1)+xx(k3)
  a2i=xx(k1+1)+xx(k3+1)
  a3r=xx(k1)-xx(k3)
  a3i=xx(k1+1)-xx(k3+1)
  xx(k0)=a0r+a2r
  xx(k0+1)=a0i+a2i
  xx(k1)=a0r-a2r
  xx(k1+1)=a0i-a2i
  xx(k2)=a1r-a3i
  xx(k2+1)=a1i+a3r
  xx(k3)=a1r+a3i
  xx(k3+1)=a1i-a3r
140 CONTINUE
GO TO 200

150  lblok2=nx
l2blok=lblok2-1
DO 160 k0=1,l2blok,2
  k1=k0+lblok2
  a1r=xx(k1)
  a1i=xx(k1+1)
  xx(k1)=xx(k0)-a1r
  xx(k1+1)=xx(k0+1)-a1i
  xx(k0)=xx(k0)+a1r
  xx(k0+1)=xx(k0+1)+a1i
160 CONTINUE


200  DO 300 m=lstart,n,3
  lblok2=nx/2**(m+1)
  l2blok=lblok2-1
  lblok1=l2blok-1
  lblok8=lblok2*8
  lblast=nx2-lblok8+1
  
  DO 210 k=4,n
    nbit(k)=0
  210 CONTINUE
  
  nw=0
  
 DO 290 offset=1,lblast,lblok8
    IF(offset == 1)GO TO 220
    arg=con1*dfloat(nw)
    w(1)=DCOS(arg)
    w(2)=DSIN(arg)
    cssqa=w(1)*w(1)
    w(3)=cssqa+cssqa-1.0D0
    w(4)=w(1)*w(2)
    w(4)=w(4)+w(4)
    w(5)=w(3)*w(1)-w(4)*w(2)
    w(6)=w(4)*w(1)+w(3)*w(2)
    cssq2a=w(3)*w(3)
    w(7)=cssq2a+cssq2a-1.0D0
    w(8)=w(4)*w(3)
    w(8)=w(8)+w(8)
    w(9)=w(7)*w(1)-w(8)*w(2)
    w(10)=w(8)*w(1)+w(7)*w(2)
    cssq3a=w(5)*w(5)
    w(11)=cssq3a+cssq3a-1.0D0
    w(12)=w(6)*w(5)
    w(12)=w(12)+w(12)
    w(13)=w(7)*w(5)-w(8)*w(6)
    w(14)=w(8)*w(5)+w(7)*w(6)
    220  lbloko=offset+lblok1
    
    DO 260 k0=offset,lbloko,2
      k1=k0+lblok2
      k2=k1+lblok2
      k3=k2+lblok2
      k4=k3+lblok2
      k5=k4+lblok2
      k6=k5+lblok2
      k7=k6+lblok2
      xk0wr=xx(k0)
      xk0wi=xx(k0+1)
      IF(offset /= 1) GO TO 240
      xk1wr=xx(k1)
      xk1wi=xx(k1+1)
      xk2wr=xx(k2)
      xk2wi=xx(k2+1)
      xk3wr=xx(k3)
      xk3wi=xx(k3+1)
      xk4wr=xx(k4)
      xk4wi=xx(k4+1)
      xk5wr=xx(k5)
      xk5wi=xx(k5+1)
      xk6wr=xx(k6)
      xk6wi=xx(k6+1)
      xk7wr=xx(k7)
      xk7wi=xx(k7+1)
      GO TO 250
      240  xk1wr=xx(k1)*w(1)-xx(k1+1)*w(2)
      xk1wi=xx(k1)*w(2)+xx(k1+1)*w(1)
      xk2wr=xx(k2)*w(3)-xx(k2+1)*w(4)
      xk2wi=xx(k2)*w(4)+xx(k2+1)*w(3)
      xk3wr=xx(k3)*w(5)-xx(k3+1)*w(6)
      xk3wi=xx(k3)*w(6)+xx(k3+1)*w(5)
      xk4wr=xx(k4)*w(7)-xx(k4+1)*w(8)
      xk4wi=xx(k4)*w(8)+xx(k4+1)*w(7)
      xk5wr=xx(k5)*w(9)-xx(k5+1)*w(10)
      xk5wi=xx(k5)*w(10)+xx(k5+1)*w(9)
      xk6wr=xx(k6)*w(11)-xx(k6+1)*w(12)
      xk6wi=xx(k6)*w(12)+xx(k6+1)*w(11)
      xk7wr=xx(k7)*w(13)-xx(k7+1)*w(14)
      xk7wi=xx(k7)*w(14)+xx(k7+1)*w(13)
      250  a0r=xk0wr+xk4wr
      a0i=xk0wi+xk4wi
      a1r=xk1wr+xk5wr
      a1i=xk1wi+xk5wi
      a2r=xk2wr+xk6wr
      a2i=xk2wi+xk6wi
      a3r=xk3wr+xk7wr
      a3i=xk3wi+xk7wi
      a4r=a0r+a2r
      a4i=a0i+a2i
      a5r=a0r-a2r
      a5i=a0i-a2i
      a6r=a1r+a3r
      a6i=a1i+a3i
      a7r=a3i-a1i
      a7i=a1r-a3r
      xx(k0)=a4r+a6r
      xx(k0+1)=a4i+a6i
      xx(k1)=a4r-a6r
      xx(k1+1)=a4i-a6i
      xx(k2)=a5r+a7r
      xx(k2+1)=a5i+a7i
      xx(k3)=a5r-a7r
      xx(k3+1)=a5i-a7i
      a0r=xk0wr-xk4wr
      a0i=xk0wi-xk4wi
      a8r=xk1wr-xk5wr
      a8i=xk1wi-xk5wi
      a1r=a8r-a8i
      a1i=a8r+a8i
      a2r=xk6wi-xk2wi
      a2i=xk2wr-xk6wr
      a8r=xk3wr-xk7wr
      a8i=xk3wi-xk7wi
      a3r=a8r-a8i
      a3i=a8r+a8i
      a4r=a0r+a2r
      a4i=a0i+a2i
      a5r=a0r-a2r
      a5i=a0i-a2i
      a6r=(a1r-a3i)/root2
      a6i=(a1i+a3r)/root2
      a7r=(a3r-a1i)/root2
      a7i=(a3i+a1r)/root2
      xx(k4)=a4r+a6r
      xx(k4+1)=a4i+a6i
      xx(k5)=a4r-a6r
      xx(k5+1)=a4i-a6i
      xx(k6)=a5r+a7r
      xx(k6+1)=a5i+a7i
      xx(k7)=a5r-a7r
      xx(k7+1)=a5i-a7i
   260 CONTINUE
    
    DO 280 k=4,n
      IF(nbit(k) /= 0)GO TO 270
      nbit(k)=1
      nw=nw+jnt(k)
      GOTO 290
      270  nbit(k)=0
      nw=nw-jnt(k)
   280 CONTINUE
    
   290 CONTINUE
   300 CONTINUE


nw=0

DO 310 k=1,n
  jnt(k)=jnt(k)+jnt(k)
  nbit(k)=0
310 CONTINUE

k=0
IF(nw <= k)GO TO 320
holdr=xx(nw+1)
holdi=xx(nw+2)
xx(nw+1)=xx(1)
xx(nw+2)=xx(2)
xx(1)=holdr
xx(2)=holdi

320  DO 340 m=1,n
  IF(nbit(m) /= 0)GO TO 330
  nbit(m)=1
  nw=nw+jnt(m)
  GO TO 350
  330  nbit(m)=0
  nw=nw-jnt(m)
340 CONTINUE

350  DO 390 k=2,nx2ls2,2
  IF(nw <= k)GO TO 360
  holdr=xx(nw+1)
  holdi=xx(nw+2)
  xx(nw+1)=xx(k+1)
  xx(nw+2)=xx(k+2)
  xx(k+1)=holdr
  xx(k+2)=holdi
  
  360  DO 380 m=1,n
    IF(nbit(m) /= 0)GO TO 370
    nbit(m)=1
    nw=nw+jnt(m)
    GO TO 390
    370  nbit(m)=0
    nw=nw-jnt(m)
  380 CONTINUE
  
390 CONTINUE

IF(signi > 0.0)GO TO 420

DO 410 k=1,nx2ls1,2
  xx(k+1)=-xx(k+1)
410 CONTINUE


420  RETURN
END SUBROUTINE dcool
