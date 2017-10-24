      program spline 

* Define program constants

      parameter (MAXPTS=100)

* Define program variables

      double precision xpt(MAXPTS), ypt(MAXPTS) 
      double precision ACOEF(MAXPTS), BCOEF(MAXPTS)
      double precision CCOEF(MAXPTS), DCOEF(MAXPTS)

* Define program functions

      double precision sfunc 
      external sfunc

*  Initilize all Array Elements to Zero
      
      call init ( MAXPTS, xpt, ypt, ACOEF, BCOEF, CCOEF, DCOEF)


*  Read in values for function

      call readdt ( N, xpt, ypt )

      call nspline ( N, xpt, ypt, ACOEF, BCOEF, CCOEF, DCOEF )

      do 87 kk=1, N
          print *, xpt(kk), ' ', ypt(kk)
87    continue

      call output ( N, MAXPTS,  xpt, ypt, ACOEF, BCOEF, CCOEF, DCOEF )
 
      end

      
      subroutine init (MAXPTS, xpt, ypt, ACOEF, BCOEF, CCOEF, DCOEF)
      double precision xpt(MAXPTS), ypt(MAXPTS) 
      double precision ACOEF(MAXPTS), BCOEF(MAXPTS)
      double precision CCOEF(MAXPTS), DCOEF(MAXPTS)
      do 5 i = 1, MAXPTS
         xpt(i) = 0.0D0
         ypt(i) = 0.0D0
         ACOEF(i) = 0.0D0
         BCOEF(i) = 0.0D0
         CCOEF(i) = 0.0D0
         DCOEF(i) = 0.0D0
5     continue
      end

      subroutine readdt (N, xpt, ypt )

      double precision xpt(*), ypt(*)
      double precision xtmp, ytmp

      open (unit=15, file='xy_pts.dat', status='old')

      N = 0 
100   read (15,*,end=200) xtmp, ytmp
         N = N + 1 
         xpt(N)=xtmp
         ypt(N)=ytmp
         goto 100
200   continue
      close(15)

      end


      subroutine nspline ( N, xpt, ypt, ACOEF, BCOEF, CCOEF, DCOEF )

      parameter (MAXPTS = 100 )

      double precision xpt(*), ypt(*) 
      double precision ACOEF(*), BCOEF(*)
      double precision CCOEF(*), DCOEF(*)
      double precision A(MAXPTS), ASUB(MAXPTS), ASUP(MAXPTS)
      double precision ASUP2(MAXPTS), B(MAXPTS)
      double precision h(MAXPTS)
      integer IPIV(MAXPTS), INFO

* Fill the diagonal, superdiagonal and subdiagonal matrix elements  
* for the natural cubic spline.  To save time, do as much work as 
* possible in only two loops (BUILD MATRIX A)  
  
      do 20 i = 1, N-1
         h(i) = xpt(i+1) - xpt(i)
         ASUB(i) = h(i) 
         ASUP(i) = h(i)
20    continue

      do 30 i = 2, N 
         A(i) = 2.0D0*(h(i-1)+h(i))  
30    continue

* Now fix all the endpoints of the arrays

      A(1) = 1.0D0
      A(N) = 1.0D0
      ASUP(1)   = 0.0D0
      ASUB(N-1) = 0.0D0

* Now build the B MATRIX

      do 40 i = 2, N-1
        B(i) = 3.0D0* ( ( ypt(i+1)-ypt(i) ) / h(i) -
     +                  ( ypt(i)-ypt(i-1) ) / h(i-1) ) 
40    continue

* Now fix all the endpoints of the array

      B(1) = 0.0D0
      B(N) = 0.0D0

* Call the appropriate LAPACK routines to solve the tridiagonal system

      CALL DGTTRF( N, ASUB, A, ASUP, ASUP2, IPIV, INFO )
      CALL DGTTRS( 'N', N, 1, ASUB, A, ASUP, ASUP2, IPIV, B, N, INFO)  
      
* Now compute the values for the various coefficient arrays

      do 50 i = 1, N-1
        ACOEF(i) = ypt(i)
        CCOEF(i) = B(i)
50    continue 

      do 60 i = N-1, 1, -1
        BCOEF(i) = ( ACOEF(i+1)-ACOEF(i) ) / h(i) -
     +           h(i) * ( CCOEF(i+1) + 2.0D0 * CCOEF(i) ) / 3.0D0 
        DCOEF(i) = ( CCOEF(i+1) - CCOEF(i) ) / (3.0D0*h(i))
60    continue

      end


      subroutine output ( N, MAXPTS, xpt, ypt, A, B, C, D) 

      double precision x, y, xpt(*), ypt(*), A(*), B(*), C(*), D(*)
      double precision delx
      double precision sfunc
      external sfunc 

      delx = xpt(N)-xpt(1) 
      delx = delx / dble(MAXPTS-1)  
      x = xpt(1)
      print *, 'N = ', N, '  MAXPTS = ', MAXPTS, ' X = ', x  
      open(unit=15, file='fit.dat', status='unknown')

      do 10 i = 1, MAXPTS
         y = sfunc(x,N,xpt,A,B,C,D)
         write (15,*) x, '  ', y
         x = x+delx
10    continue
      close(15)
      end       
       
     

      double precision function sfunc (x,N,xpt,A,B,C,D)

      double precision x, xpt(*), A(*), B(*), C(*), D(*)
      integer loc

      external loc

*     find location of x in arrays using bisection

      i = loc(xpt,N,x)
      sfunc = A(i) + B(i) * (x - xpt(i)) + C(i) * ( x - xpt(i))**2
     +               + D(i) * ( x - xpt(i))**3
      return
      end



      integer function loc(xpt,n,x)
      integer n
      double precision x,xpt(n)
      integer spindex
      logical found
      found = .false.
      spindex = -1
      do i=1,n-1
        if ( (.not. found) .and. 
     +       ((x .ge. xpt(i)) .and. (x .lt. xpt(i+1))) ) then
             spindex = i
             found = .true.
        endif
        print *, s, xpt(i), xpt(i+1), spindex
      enddo
      if ( .not. found ) spindex = n
      loc=spindex
      return
      end 








