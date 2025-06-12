!**********************************************************************************************************************************
!
!                                                          L I N R E G
!
!  Program:      LINREG
!
!  Programmer:   Dr. David G. Simpson
!                Department of Physical Science
!                Prince George's Community College
!                Largo, Maryland  20774
!
!  Date:         January 21, 2002
!
!  Language:     Fortran-90
!
!  Description:  This program performs a linear regression analysis for a set of data given as (x,y) pairs.  The output from
!                the program is the slope and y-intercept of the least-squares best fit straight line through the data points.
!
!**********************************************************************************************************************************


!**********************************************************************************************************************************
!  Main program
!**********************************************************************************************************************************

subroutine linreg_sub(x,y,n,b,m,r)
!program linreg

   implicit none                                                                    ! no default data types

   integer, parameter  :: dbl = kind (0.0d0)                                        ! define kind for double precision

   real(dbl), intent(out)           ::  b                                                        ! y-intercept of least-squares best fit line
   real(dbl), intent(out)           ::  m                                                        ! slope of least-squares best fit line
   real( kind=8 ), intent(out)           ::  r                                                       ! squared correlation coefficient

   integer, intent(in)           ::  n                                              !INPUT     ! number of data points
   real, intent(in)           ::  x(n)                                           !INPUT     ! input x data
   real, intent(in)           ::  y(n)                                           !INPUT     ! input y data

   integer :: i
   character (len=80)  ::  str                                                      ! input string
   real(dbl)           :: sumx                                            ! sum of x
   real(dbl)           :: sumx2                                           ! sum of x**2
   real(dbl)           :: sumxy                                           ! sum of x * y
   real(dbl)           :: sumy                                            ! sum of y
   real(dbl)           :: sumy2                                           ! sum of y**2

!   write (unit=*, fmt="(a)") " LINREG - Perform linear regression"                  ! print introductory message
!   write (unit=*, fmt="(a/)") "   (Enter END to stop data entry and compute"//  &
!                              " linear regression.)"

!Initialisation
b = 0.0
m = 0.0
r = 0.0
sumx  = 0.0
sumx2 = 0.0
sumxy = 0.0
sumy  = 0.0
sumy2 = 0.0

   do 100 i=1,n,1                                                                              ! loop for all data points

!      write (unit=*, fmt="(a)", advance="no") " Enter x y:  "                       ! prompt to enter data
!      read (unit=*, fmt="(a)") str                                                  ! read x and y into string
!      if (str == "end" .or. str == "END") exit                                      ! if no more data, then exit loop
!      read (unit=str, fmt=*) x, y                                                   ! else read x and y from string
!
!      n = n + 1.0d0                                                                 ! increment number of data points by 1
      sumx  = sumx + x(i)                                                              ! compute sum of x
      sumx2 = sumx2 + x(i) * x(i)                                                         ! compute sum of x**2
      sumxy = sumxy + x(i) * y(i)                                                         ! compute sum of x * y
      sumy  = sumy + y(i)                                                              ! compute sum of y
      sumy2 = sumy2 + y(i) * y(i)                                                         ! compute sum of y**2

!      print*, 'i',i,'x',x(i),'y',y(i), 'sumx2',sumx2,'sumy2',sumy2,'sumxy',sumxy


100 continue
!   end do

   m = (n * sumxy  -  sumx * sumy) / (n * sumx2 - sumx**2)                          ! compute slope
   b = (sumy * sumx2  -  sumx * sumxy) / (n * sumx2  -  sumx**2)                    ! compute y-intercept
   r = (sumxy - sumx * sumy / n) /                                     &            ! compute correlation coefficient
                     sqrt((sumx2 - sumx**2/n) * (sumy2 - sumy**2/n))

!   write (unit=*, fmt="(/a,es15.6)") " Slope        m = ", m                        ! print results
!   write (unit=*, fmt="(a, es15.6)") " y-intercept  b = ", b
!   write (unit=*, fmt="(a, es15.6)") " Correlation  r = ", r

return
!end 
end subroutine linreg_sub
