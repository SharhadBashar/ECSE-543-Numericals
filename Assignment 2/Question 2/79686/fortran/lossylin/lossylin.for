C*************************************************************
C***********                                       ***********
C*********** One-dimensional demonstration program ***********
C***********                                       ***********
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C
C      The subroutines  that make up this program  communicate
C      via named common blocks.  The variables in commons are:
C
C      Problem definition and solution
C      common /problm/
C        nodes   =  number of nodes in problem
C        resist  =  resistance per unit length
C        conduc  =  conductance per unit length
C        x       =  nodal coordinates
C        voltag  =  nodal voltage values
C
C      Global matrices and right-hand side
C      common /matrix/
C        stcon   =  global coefficient matrix (connected)
C        rthdsd  =  right-hand side of system of equations
C        c       =  connection matrix for whole problem
C        stdis   =  global coefficient matrix (disconnected)
C
C      Predefined problem size parameter (array dimensions):
C        MAXNOD  =  maximum number of nodes in problem
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 30)
      common /problm/ nodes, resist, conduc,
     *                x(MAXNOD), voltag(MAXNOD)
      common /matrix/ stcon(MAXNOD,MAXNOD), rthdsd(MAXNOD),
     *      c(MAXNOD,2*MAXNOD-2), stdis(2*MAXNOD-2,2*MAXNOD-2)
C=============================================================
C
C     Initialize matrices to zeros.
      call matini(stcon, MAXNOD, MAXNOD)
      call matini(stdis, 2*MAXNOD-2, 2*MAXNOD-2)
      call vecini(rthdsd, MAXNOD)
C
C     Problem definition:
C            Sending-end voltage, resistivity, conductivity
C            Node number, x-coordinate (1 or more times)
      call prblem
C
C     Create connection matrix
      call makec
C
C     Construct the disjoint coefficient matrix
      call disjnt
C
C     Construct connected coefficient matrix
      call connct
C
C     Construct right-hand side vector
      call rtsd
C
C     Solve the assembled finite element equations
      call eqsolv(stcon, voltag, rthdsd, nodes-1, MAXNOD)
C
C     Write the resulting voltage values to output file
      call output
C
      stop
      end
C
C*************************************************************
C
      Subroutine prblem
C
C*************************************************************
C
C     Reads in problem details from input file.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 30)
      common /problm/ nodes, resist, conduc,
     *                x(MAXNOD), voltag(MAXNOD)
      common /matrix/ stcon(MAXNOD,MAXNOD), rthdsd(MAXNOD),
     *      c(MAXNOD,2*MAXNOD-2), stdis(2*MAXNOD-2,2*MAXNOD-2)
C=============================================================
C
C     Set up initial values
      nodes = 1
      do 10 i = 1,MAXNOD
   10    x(i) = 0.
C
C     Read problem data
C        (1) sending voltage, resistivity, conductivity
      read (*, *) sendgv, resist, conduc
C        (2) one or more nodes:  node number, x-coordinate
   30 read (*, *, end=60) node, xnode
      if (node .le. 0 .or. nodes+node .gt. MAXNOD) then
	call errexc('PRBLEM', 1)
      else
	xleft = x(nodes)
	nleft = nodes
	xnode = xleft + xnode
	nodes = nodes + node
	segms = nodes - nleft
	do 40 node=nleft,nodes
	  x(node) = (nodes-node)*xleft - (nleft-node)*xnode
   40     x(node) = x(node) / segms
      endif
      go to 30
C
C     Data successfully read.  Set up problem values.
   60 if (nodes .le. 1) then
	call errexc('PRBLEM', 2)
      endif
      if (sendgv .eq. 0.) sendgv = 1.
      if (resist .eq. 0.) resist = 1.
      if (conduc .eq. 0.) conduc = 1.
      do 70 i = 1,nodes
   70   voltag(i) = 0.
      voltag(nodes) = sendgv
C
      return
      end
C
C*************************************************************
C
      Subroutine makec
C
C*************************************************************
C
C     Establishes connection matrix c.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 30)
      common /problm/ nodes, resist, conduc,
     *                x(MAXNOD), voltag(MAXNOD)
      common /matrix/ stcon(MAXNOD,MAXNOD), rthdsd(MAXNOD),
     *      c(MAXNOD,2*MAXNOD-2), stdis(2*MAXNOD-2,2*MAXNOD-2)
C=============================================================
C
C     Create connection matrix:
      do 60 i = 1,nodes
	do 50 j = 1,2*nodes-2
	  if (j/2+1 .eq. i) then
	    c(i,j) = 1.
	  else
	    c(i,j) = 0.
	  endif
   50     continue
   60   continue
C
      return
      end
C
C*************************************************************
C
      Subroutine disjnt
C
C*************************************************************
C
C     Constructs the disjoint coefficient matrix
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 30)
      common /problm/ nodes, resist, conduc,
     *                x(MAXNOD), voltag(MAXNOD)
      common /matrix/ stcon(MAXNOD,MAXNOD), rthdsd(MAXNOD),
     *      c(MAXNOD,2*MAXNOD-2), stdis(2*MAXNOD-2,2*MAXNOD-2)
C=============================================================
C
      do 70 i = 1,2*nodes-3,2
	j = i + 1
C
C        Determine element length, quit if nonpositive
	ellng = x(j/2+1) - x(j/2)
	if (ellng .le. 0.) then
	  call errexc('DISJNT', i)
	else
C           Fit element s and t into disjoint global matrix
	  stdis(i,i) = + 1./(resist*ellng) + conduc*ellng/3.
	  stdis(i,j) = - 1./(resist*ellng) + conduc*ellng/6.
	  stdis(j,i) = - 1./(resist*ellng) + conduc*ellng/6.
	  stdis(j,j) = + 1./(resist*ellng) + conduc*ellng/3.
	endif
   70   continue
C
      return
      end
C
C*************************************************************
C
      Subroutine connct
C
C*************************************************************
C
C     Connection transformation:     stcon  =  c' stdis c
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 30)
      common /problm/ nodes, resist, conduc,
     *                x(MAXNOD), voltag(MAXNOD)
      common /matrix/ stcon(MAXNOD,MAXNOD), rthdsd(MAXNOD),
     *      c(MAXNOD,2*MAXNOD-2), stdis(2*MAXNOD-2,2*MAXNOD-2)
C=============================================================
C
C     Set connected coefficient matrix to zero:
      do 50 k = 1,nodes
	do 40 l = 1,nodes
	  sum = 0.
	  do 30 j = 1,2*nodes-2
	    do 20 i = 1,2*nodes-2
	      sum = sum + c(k,i) * stdis(i,j) * c(l,j)
   20         continue
   30       continue
	  stcon(k,l) = sum
   40     continue
   50   continue
C
      return
      end
C
C*************************************************************
C
      Subroutine rtsd
C
C*************************************************************
C
C     Transposes sending-end voltage to right-hand side
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 30)
      common /problm/ nodes, resist, conduc,
     *                x(MAXNOD), voltag(MAXNOD)
      common /matrix/ stcon(MAXNOD,MAXNOD), rthdsd(MAXNOD),
     *      c(MAXNOD,2*MAXNOD-2), stdis(2*MAXNOD-2,2*MAXNOD-2)
C=============================================================
C
      do 10 i = 1,nodes-1
	rthdsd(i) = -stcon(i,nodes) * voltag(nodes)
   10   continue
      return
      end
C*************************************************************
C
      Subroutine output
C
C*************************************************************
C
C     Outputs problem and results to standard output stream,
C     accompanied by exact analytic solution.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 30)
      common /problm/ nodes, resist, conduc,
     *                x(MAXNOD), voltag(MAXNOD)
      common /matrix/ stcon(MAXNOD,MAXNOD), rthdsd(MAXNOD),
     *      c(MAXNOD,2*MAXNOD-2), stdis(2*MAXNOD-2,2*MAXNOD-2)
C=============================================================
C
      p = sqrt(resist * conduc)
      dn = exp(p*x(nodes)) + exp(-p*x(nodes))
C
      nbits = 100
      do 10 i = 1,nbits+1
	z = (i-1) * x(nodes)/nbits
	vap = vapx(z)
	vex = voltag(nodes) * (exp(p*z) + exp(-p*z))/dn
	err = (vap - vex)/vex * 100.
	write (*, 1000, err=890) i, z, vap, vex, err
   10   continue
 1000 format (1x, i3, 4f12.6)
      return
  890 call errexc('OUTPUT', 1)
      end
C
C*************************************************************
C
      Function vapx(z)
C
C*************************************************************
C
C     Returns interpolated approximate solution vapx at x = z
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 30)
      common /problm/ nodes, resist, conduc,
     *                x(MAXNOD), voltag(MAXNOD)
      common /matrix/ stcon(MAXNOD,MAXNOD), rthdsd(MAXNOD),
     *      c(MAXNOD,2*MAXNOD-2), stdis(2*MAXNOD-2,2*MAXNOD-2)
C=============================================================
C
C     Determine in which interval z lies
      int = 0
      do 10 i = 1,nodes-1
	 if (x(i) .le. z .and. x(i+1) .ge. z) int = i
   10    continue
      if (int .eq. 0) then
	call errexc('VAPX', 1)
      else
C         Interpolate within interval to find value
	vapx = (voltag(int+1) * (z - x(int)) +
     *     voltag(int) * (x(int+1) - z)) / (x(int+1) - x(int))
      endif
C
      return
      end
