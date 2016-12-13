C*************************************************************
C********                                            *********
C******** Strip and ground plane:  integral elements *********
C********                                            *********
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C
C      The subroutines  that make up this program  communicate
C      via named common blocks.  The variables in commons are:
C
C      Problem definition and solution
C      common /problm/
C        space   =  spacing, strip to ground plane
C        nelmt   =  number of elements in half-strip
C        x       =  nodal coordinates
C        charg   =  charge density on substrip
C
C      Global matrices and right-hand side
C      common /matrix/
C        coeffs   =  global coefficient matrix (connected)
C        rthdsd  =  right-hand side of system of equations
C
C      Predefined problem parameters (array dimensions):
C        MAXELM  =  maximum number of elements in problem
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXELM = 100, PI = 3.141596)
      common /problm/ space, nelmt, x(MAXELM+1), charg(MAXELM)
      common /matrix/ coeffs(MAXELM,MAXELM), rthdsd(MAXELM)
C=============================================================
C
C     Initialize matrices to zeros.
      call matini(coeffs, MAXELM, MAXELM)
      call vecini(charg, MAXELM)
C
C     Problem definition:
C               spacing between strip and ground plane
C               number of elements, x-value (1 or more times)
      call inprob
C
C     Create coefficient matrix
      call mkcoef
C
C     Construct right-hand side vector
      do 20 j = 1,nelmt
   20   rthdsd(j) = 1.
C
C     Solve the assembled finite element equations
      call eqsolv(coeffs, charg, rthdsd, nelmt, MAXELM)
C
C     Write the resulting charge densities to output file
      call output
C
      stop
      end
C
C*************************************************************
C
      Subroutine inprob
C
C*************************************************************
C
C     Reads in problem details from input file.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXELM = 100, PI = 3.141596)
      common /problm/ space, nelmt, x(MAXELM+1), charg(MAXELM)
      common /matrix/ coeffs(MAXELM,MAXELM), rthdsd(MAXELM)
C=============================================================
C
C     Set up initial values
      nelmt = 1
      call vecini(x, MAXELM+1)
C
C     Read problem data
C        (1) spacing from ground plane to strip
      read (*, *) space
C        (2) number of elements, x at end of element group
C            (may be repeated any number of times)
   30 read (*, *, end=60) nelm, xelmt
      if (nelm .le. 0 .or. nelmt+nelm .gt. MAXELM) then
	call errexc('INPROB', 1000 + nelm)
      else
	xleft = x(nelmt)
	nleft = nelmt
	xelmt = xleft + xelmt
	nelmt = nelmt + nelm
	segms = nelmt - nleft
	do 40 nelm=nleft,nelmt
	  x(nelm) = (nelmt-nelm)*xleft - (nleft-nelm)*xelmt
   40     x(nelm) = x(nelm) / segms
      endif
      go to 30
C
C     Data successfully read.  Set up problem values.
   60 nelmt = nelmt - 1
      if (nelmt .le. 1) call errexc('INPROB', 2)
C
      return
      end
C
C*************************************************************
C
      Subroutine mkcoef
C
C*************************************************************
C
C     Establishes coefficient matrix coeffs.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXELM = 100, PI = 3.141596)
      common /problm/ space, nelmt, x(MAXELM+1), charg(MAXELM)
      common /matrix/ coeffs(MAXELM,MAXELM), rthdsd(MAXELM)
C=============================================================
      do 80 i = 1,nelmt
C
C       Determine element i width and midpoint, quit if < 0
	widti = x(i+1) - x(i)
	xmidi = (x(i+1) + x(i))/2.
	if (widti .le. 0.) then
	  call errexc('MKCOEF', i)
	else
C         Compute entries of coeffs from Green's functions
	  do 70 j = 1,i
	    xmidj = (x(j+1) + x(j))/2.
	    dist1 = abs(xmidi - xmidj)
	    if (i .eq. j) dist1 = exp(alog(widti) - 1.5)
	    dist2 = xmidi + xmidj
	    dist3 = sqrt((xmidi + xmidj)**2 + (2*space)**2)
	    dist4 = sqrt((xmidi - xmidj)**2 + (2*space)**2)
	    coeffs(i,j) = alog((dist1*dist2)/(dist3*dist4))
     *                         / (-2. * PI)
	    coeffs(j,i) = coeffs(i,j)
   70       continue
	endif
   80   continue
C
      return
      end
C
C*************************************************************
C
      Subroutine output
C
C*************************************************************
C
C     Outputs problem and results to default output stream.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXELM = 100, PI = 3.141596)
      common /problm/ space, nelmt, x(MAXELM+1), charg(MAXELM)
      common /matrix/ coeffs(MAXELM,MAXELM), rthdsd(MAXELM)
C=============================================================
C
      totchg = 0.
      do 10 i = 1,nelmt
	width = x(i+1) - x(i)
	denst = charg(i) / width
	totchg = totchg + charg(i)
	write (*, 1000, err=890) i, x(i), denst
	write (*, 1000, err=890) i+1, x(i+1), denst
   10   continue
 1000 format (1x, i3, 2f12.6)
      write (*, 1000, err=890) nelmt, totchg
      return
  890 call errexc('OUTPUT', 1)
      end
