C*************************************************************
C
      Subroutine eqsolv(s, x, y, n, maxn)
C
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C     Solves the system of equations with square symmetric co-
C     efficient matrix s,  in n variables but dimensions maxn,
C                          s * x  =  y
C     Procedure:   Cholesky decomposition, forward elimination
C     and back substitution.  s and y are not preserved.
C=============================================================
      dimension s(maxn,maxn), x(maxn), y(maxn)
C
C     Cholesky decomposition replaces s by triangular factors
      call chlsky(s, s, n, maxn)
C         then forward eliminate and back substitute.
      call fwdelm(s, y, y, n, maxn)
      call backsb(s, y, x, n, maxn)
      return
      end
C
C*************************************************************
C
      Subroutine chlsky(a, fl, n, maxn)
C
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
      dimension a(maxn,maxn), fl(maxn,maxn)
C
      do 70 j = 1,n
	fl(j,j) = a(j,j)
	do 20 k = 1,j-1
	  fl(j,j) = fl(j,j) - fl(j,k)*fl(j,k)
   20     continue
	if (fl(j,j) .gt. 0.) then
	  fl(j,j) = sqrt(fl(j,j))
	  do 50 i = j+1,n
	    fl(i,j) = a(i,j)
	    do 40 k = 1,j-1
	      fl(i,j) = fl(i,j) - fl(i,k)*fl(j,k)
   40         continue
	    fl(i,j) = fl(i,j) / fl(j,j)
	    fl(j,i) = fl(i,j)
   50       continue
	else
	  call errexc('CHLSKY', j)
	endif
   70   continue
      return
      end
C
C*************************************************************
C
      Subroutine fwdelm(fl, b, y, n, maxn)
C
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C     Forward elimination of triangular system  (fl)*y = b
C=============================================================
      dimension fl(maxn,maxn), b(maxn), y(maxn)
C
      y(1) = b(1) / fl(1,1)
      do 60 i = 2,n
	y(i) = b(i)
	do 40 j = 1,i-1
	  y(i) = y(i) - fl(i,j)*y(j)
   40     continue
	y(i) = y(i) / fl(i,i)
   60   continue
      return
      end
C
C*************************************************************
C
      Subroutine backsb(u, b, x, n, maxn)
C
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C     Back substitution of triangular system  u*x = b
C=============================================================
      dimension u(maxn,maxn), b(maxn), x(maxn)
C
      x(n) = b(n) / u(n,n)
      do 60 i = n-1,1,-1
	x(i) = b(i)
	do 40 j = i+1,n
	  x(i) = x(i) - u(i,j)*x(j)
   40     continue
	x(i) = x(i) / u(i,i)
   60   continue
      return
      end
C
C*************************************************************
C
      Subroutine vecini(a, nr)
C
C*************************************************************
C
C     Initializes vector a, of nr entries, to all zeros.
C=============================================================
      dimension a(nr)
C=============================================================
C
      do 10 i = 1,nr
	a(i) = 0.
   10   continue
      return
      end
C
C*************************************************************
C
      Subroutine matini(a, nr, nc)
C
C*************************************************************
C
C     Initializes  matrix a, of nr rows and nc columns, to all
C     zero entries.
C=============================================================
      dimension a(nr, nc)
C=============================================================
C
      do 20 i = 1,nr
	do 10 j = 1,nc
	  a(i,j) = 0.
   10     continue
   20   continue
      return
      end
C
C*************************************************************
C
      Function locate(i,j)
C
C*************************************************************
C
C     Returns linear storage-mode  location index of the (i,j)
C     element of  a symmetric matrix  whose lower  triangle is
C     stored by rows in a linear array.
C
      if (j .lt. i) then
	locate = (i*(i-1))/2+j
      else
	locate = (j*(j-1))/2+i
      endif
      return
      end
C
C*************************************************************
C
      Subroutine errexc(prog, ierr)
C
C*************************************************************
C
C     Exception handler for serious errors: prints the program
C     name and error number, then halts execution.
C
C     Execution is not halted  if the logical variable "debug"
C     is true; this is often useful for program testing.
C
      Character*(*) prog
      Logical debug
      debug = .true.
      if (debug) then
	 write (*,100) prog, ierr
	 return
      else
	 write (*,101) prog, ierr
	 stop 911
      endif
  100 format(1x/' Error  in  routine ', a6, ', error ',
     *       'number ', i5)
  101 format(1x/' Stopped in routine ', a6, ', error ',
     *       'number ', i5)
      end
