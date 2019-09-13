      program linterp

* Define program constants

      parameter (MAXPTS=1000)
      parameter (MAXVAL=1000)

* Define program variables

      double precision xpt(0:MAXPTS), ypt(0:MAXPTS)
      double precision x(0:MAXPTS)
      double precision xtmp, ytmp, xmin, xmax, delx


* Define program functions

      double precision interp
      external interp

*  Read in values for function

      open (unit=15, file='xy_pts.dat', status='old')

      numpts = -1 
100   read (15,*,end=200) xtmp, ytmp
         numpts=numpts+1
         xpt(numpts)=xtmp
         ypt(numpts)=ytmp
* If you need to see the input data, ucomment the following line
*         print *, numpts, '  ', xpt(numpts), '  ', ypt(numpts)
         goto 100
200   continue
      close(15)

*  Find maximum and minimum in domain

      xmin=xpt(0)
      xmax=xpt(0)
      do 10 i = 1, numpts
        xmin = min(xmin,xpt(i))
        xmax = max(xmax,xpt(i))
10    continue

* Set values for domain of interpolating polynomial

      delx = (xmax-xmin)/dble(MAXVAL-1) 
      x(1)=xmin
      do 20 i = 2, MAXVAL
         x(i) = x(i-1) + delx
20    continue

* Nom compute values of interpolating polynomial

      do 30 i = 1, MAXVAL
         print *, x(i),  interp(numpts, xpt, ypt, x(i) )
30    continue

      end


      double precision function interp ( N, xpt, ypt, x )

      double precision xpt(0:*), ypt(0:*), x
      double precision sum, prod 
      
      sum = 0.0D0 
      do 10 k = 0, N
          prod =1.0D0
          do 20 i = 0, N
             if (i .ne. k) prod = prod * (x-xpt(i))/(xpt(k)-xpt(i))
20        continue
          sum = sum + ypt(k) * prod
10    continue

      interp = sum
      return
      end  


