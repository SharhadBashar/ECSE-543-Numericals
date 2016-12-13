C*************************************************************
C************                                     ************
C************      Saturable material program     ************
C************                                     ************
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C
C      Reads a saturable-material  problem from standard input
C      and writes the solution to the standard output stream.
C
C      The subroutines  that make up this program  communicate
C      via named common blocks.  The variables in commons are:
C
C      Global parameter values needed by all program units
C      common /problm/
C        nodes   =  number of nodes used in problem
C        nelmts  =  number of elements in model
C        x, y    =  nodal coordinates
C        constr  =  logical, .true. for fixed potentials
C        potent  =  nodal potential array
C        step    =  Newton step array
C        nvtx    =  list of nodes for each element
C        source  =  source density in each element
C        materl  =  material code for each element
C
C      Global matrix and right-hand side
C      common /matrix/
C        p       =  global P-matrix for whole problem
C        resid   =  residual vector (right-hand side)
C        step    =  Newton step computed from resid
C        b2max   =  largest B-squared in any element
C
C      Temporary working arrays and variables
C      common /workng/
C        pel     =  P for one element (working array)
C        sel     =  S for one element (working array)
C        tel     =  T for one element (working array)
C
C      Predefined problem size parameters (array dimensions):
C        MAXNOD  =  maximum number of nodes in problem
C        MAXELM  =  maximum number of elements in problem
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75)
      parameter (TINY = 1.e-35, HUGE = 1.e+35)
      logical constr
      common /problm/ nodes, nelmts,
     *  x(MAXNOD), y(MAXNOD), constr(MAXNOD), potent(MAXNOD),
     *  nvtx(3,MAXELM), source(MAXELM), materl(MAXELM)
      common /matrix/ p(MAXNOD,MAXNOD), resid(MAXNOD),
     *                step(MAXNOD), b2max
      common /workng/ pel(3,3), sel(3,3), tel(3,3), e(3)
C=============================================================
C
      parameter (TOLER = 1.e-5)
C
C     Fetch input data from input file
      call meshin
C
C     Newton iteration: set up to start
      Write (*,500)
  500 Format (1x/ 16x, 'Newton iteration' // 1x, 'newt', 4x,
     *  'B max', 3x, 'dampg', 3x, 'stepnorm', 4x, 'potnorm',
     *  4x,' convgc')
      convgc = 1.
      newt = 0
C
C        Do up to 12 steps while 'convgc' exceeds tolerance
    1 if (newt .lt. 12  .and.  convgc .gt. TOLER) then
	newt = newt + 1
	b2max = 0.
C       Set global p-matrix and residual to all zeros.
	call matini(p, MAXNOD, MAXNOD)
	call vecini(resid, MAXNOD)
C
C       Assemble global matrix, element by element.
	do 10 i = 1,nelmts
C           Construct element P, S and T matrices
	  call elmatr(i)
C           Embed matrices in global P; augment right side:
	  call elembd(i)
   10     continue
C
C          Solve the assembled finite element equations
	call eqsolv(p, step, resid, nodes, MAXNOD)
C          Fix up Newton step and check convergence.
C          First find vector norms
	  ptnorm = TINY
	  stnorm = TINY
	  do 20 i = 1,nodes
	    ptnorm = max(ptnorm, abs(potent(i)), abs(step(i)))
	    stnorm = max(stnorm, abs(step(i)))
   20       continue
	  convgc = stnorm/ptnorm
	  damper = 1. - 0.5 * min(1.0, convgc**2)
	  do 30 i = 1,nodes
	    potent(i) = potent(i) + damper * step(i)
   30       continue
	write (*,510) newt, sqrt(b2max), damper, stnorm,
     *                ptnorm, convgc
	go to 1
      endif
  510   Format (1x, i3, 2x, 0p2f8.3, 1p5e11.2)
C
C     Print out the resulting potential values
      call output
C
      stop
      end
C
C*************************************************************
C
      Subroutine elmatr(ie)
C
C*************************************************************
C
C     Constructs element matrices  p and t for a single first-
C     order triangular finite element.  ie = element number.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75)
      parameter (TINY = 1.e-35, HUGE = 1.e+35)
      logical constr
      common /problm/ nodes, nelmts,
     *  x(MAXNOD), y(MAXNOD), constr(MAXNOD), potent(MAXNOD),
     *  nvtx(3,MAXELM), source(MAXELM), materl(MAXELM)
      common /matrix/ p(MAXNOD,MAXNOD), resid(MAXNOD),
     *                step(MAXNOD), b2max
      common /workng/ pel(3,3), sel(3,3), tel(3,3), e(3)
C=============================================================
C
C     Set up indices for triangle
      i = nvtx(1,ie)
      j = nvtx(2,ie)
      k = nvtx(3,ie)
C
C     Compute element T-matrix
      area = abs((x(j) - x(i)) * (y(k) - y(i)) -
     1           (x(k) - x(i)) * (y(j) - y(i))) / 2.
      do 20 l = 1,3
	do 10 m = 1,3
	  tel(l,m) = area / 12.
   10     continue
	tel(l,l) = 2. * tel(l,l)
   20   continue
C
C     Compute element geometric P-matrix
      i1 = 1
      i2 = 2
      i3 = 3
      do 30 l = 1,3
	do 30 m = 1,3
	  sel(l,m) = 0.
   30     continue
      do 50 nvrtex = 1,3
	ctng = ((x(j) - x(i)) * (x(k) - x(i)) +
     1         (y(j) - y(i)) * (y(k) - y(i))) / (2. * area)
	ctng2 = ctng / 2.
	sel(i2,i2) = sel(i2,i2) + ctng2
	sel(i2,i3) = sel(i2,i3) - ctng2
	sel(i3,i2) = sel(i3,i2) - ctng2
	sel(i3,i3) = sel(i3,i3) + ctng2
C           Permute row and column indices once
	i4 = i1
	i1 = i2
	i2 = i3
	i3 = i4
	l = i
	i = j
	j = k
	k = l
   50   continue
C
C     Find squared flux density in element
      bsq = 0.
      do 70 i = 1,3
	e(i) = 0.
	do 60 j = 1,3
	  e(i) = e(i) + sel(i,j)*potent(nvtx(j,ie))
   60     continue
	bsq = bsq + e(i) * potent(nvtx(i,ie))/area
   70   continue
      b2max = max(bsq, b2max)
C
C     Find reluctivity and derivative for material not 0
      if (materl(ie) .ne. 0) then
	call reluc(bsq, rnu, dnu)
      else
	rnu = 1.
	dnu = 0.
      endif
C
C     Create Jacobian P, and S matrix for element
      do 90 i = 1,3
	do 80 j = 1,3
	  sel(i,j) = rnu * sel(i,j)
	  pel(i,j) = sel(i,j) + 2.*dnu*e(i)*e(j)/area
   80     continue
   90   continue
C
      return
      end
C
C*************************************************************
C
      Subroutine reluc(b2, rval, rder)
C
C*************************************************************
C
C     Find the  relative reluctivity and  its derivative  with
C     respect to b2, the squared flux density. The material is
C     a rather ordinary sheet steel.
C
      Dimension b0(11), r0(11), r2(11)
      Data  N, b0, r0, r2       /     11, 0.0000000E+00,
     *      0.4000000E-01, 0.1600000E+00, 0.3600000E+00,
     *      0.6400000E+00, 0.1000000E+01, 0.1440000E+01,
     *      0.1690000E+01, 0.1960000E+01, 0.2250000E+01,
     *      0.2560000E+01, 0.1090000E-03, 0.1090000E-03,
     *      0.1090000E-03, 0.1100000E-03, 0.1110000E-03,
     *      0.1120000E-03, 0.1200000E-03, 0.1430000E-03,
     *      0.1960000E-03, 0.4170000E-03, 0.1250000E-02,
     *      0.0000000E+00, 0.0000000E+00, 0.0000000E+00,
     *      0.0000000E+00, 0.4393732E-05, 0.0000000E+00,
     *      0.2483705E-03, 0.4209161E-03, 0.4664167E-03,
     *      0.9512425E-02, 0.0000000E+00/
C
      if (b2 .le. b0(N)) then
	kl = 1
	kr = N
    1   if (kr-kl .gt. 1) then
	  k = (kr+kl)/2
	  if (b0(k) .gt. b2) then
	    kr = k
	  else
	    kl = k
	  endif
	  go to 1
	endif
	dx = b0(kr)-b0(kl)
	du = (b0(kr)-b2)/dx
	dl = (b2-b0(kl))/dx
	du2 = du**2
	dl2 = dl**2
	rval = du*r0(kl)+dl*r0(kr) + ((du2-1.)*du*r2(kl)
     *         + (dl2-1.)*dl*r2(kr))*dx**2/6.
	rder = -(r0(kl)-r0(kr))/dx - ((3.*du2-1.)*r2(kl)
     *         - (3.*dl2-1.)*r2(kr))*dx/6.
      else
	dx = b0(N)-b0(N-1)
	rder = -(r0(N-1)-r0(N))/dx+(+r2(N-1)+2.*r2(N))*dx/6.
	rval = r0(N)+(b2-b0(N))*rder
      endif
      return
      end
C
C*************************************************************
C
      Subroutine elembd(ie)
C
C*************************************************************
C
C     Embeds single-element p and t matrices currently in sel
C     and tel (in general common block) in the global  matrix
C     p.  Argument ie is the element number.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75)
      parameter (TINY = 1.e-35, HUGE = 1.e+35)
      logical constr
      common /problm/ nodes, nelmts,
     *  x(MAXNOD), y(MAXNOD), constr(MAXNOD), potent(MAXNOD),
     *  nvtx(3,MAXELM), source(MAXELM), materl(MAXELM)
      common /matrix/ p(MAXNOD,MAXNOD), resid(MAXNOD),
     *                step(MAXNOD), b2max
      common /workng/ pel(3,3), sel(3,3), tel(3,3), e(3)
C=============================================================
C
C     Run through element p and t matrices (pel, sel and tel),
C     augmenting the global p and the right-hand side as
C     appropriate.
C
      do 30 i = 1,3
	irow = nvtx(i,ie)
C
C          Does row correspond to a fixed potential?
	if (constr(irow)) then
	  p(irow,irow) = 1.
	  resid(irow) = 0.
	else
C          No, potential is free to vary.  Do all 3 columns.
	  do 20 j = 1,3
	    icol = nvtx(j,ie)
C
C          Column corresponds to fixed potential? Augment
C              residual only (otherwise also augment matrix)
	    if (constr(icol)) then
	      resid(irow) = resid(irow) + tel(i,j) *
     *                    source(ie) - sel(i,j) * potent(icol)
	    else
	      p(irow,icol) = p(irow,icol) + pel(i,j)
	      resid(irow) = resid(irow) + tel(i,j) *
     *                    source(ie) - sel(i,j) * potent(icol)
	    endif
   20       continue
	endif
   30   continue
C
C     All done -- return to calling program.
      return
      end
C
C*************************************************************
C
      Subroutine meshin
C
C*************************************************************
C
C     Reads input data file in three parts:   nodes, elements,
C     fixed potentials.   Each division  is followed by a line
C     containing nothing but the / character,  which serves as
C     a terminator.
C
C     nodes:              node number, x, y.
C     elements:           node numbers, source.
C     potentials:         node number, fixed value.
C
C     All data are echoed as read but little or no checking is
C     done for validity.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75)
      parameter (TINY = 1.e-35, HUGE = 1.e+35)
      logical constr
      common /problm/ nodes, nelmts,
     *  x(MAXNOD), y(MAXNOD), constr(MAXNOD), potent(MAXNOD),
     *  nvtx(3,MAXELM), source(MAXELM), materl(MAXELM)
      common /matrix/ p(MAXNOD,MAXNOD), resid(MAXNOD),
     *                step(MAXNOD), b2max
      common /workng/ pel(3,3), sel(3,3), tel(3,3), e(3)
C=============================================================
C
      Dimension nold(3), nnew(3)
C
C     Read in the node list and echo input lines.
C        Start by printing a heading for the node list.
      write (*, 1140)
 1140 format (1x // 8x, 'Input node list' // 3x, 'n', 8x, 'x',
     1        11x, 'y' / 1x)
C
C        Read and echo nodes
      nodes = 0
      xold = HUGE
      yold = HUGE
   20 read (*,*, end=911) xnew, ynew
      if (xnew .ne. xold  .or.  ynew .ne. yold) then
	nodes = nodes + 1
	x(nodes) = xnew
	y(nodes) = ynew
	xold = xnew
	yold = ynew
	write (*, 1105) nodes, x(nodes), y(nodes)
	go to 20
      endif
 1105 format (1x, i3, 2(2x, f10.5))
C
C     Read in the element list and echo all input as received.
C          Print heading to start.
      Write (*, 1160)
 1160 format (1x // 6x, 'Input element list' // 3x, 'i', 5x,
     1        'j', 5x, 'k', 6x, 'Source', 3x, 'Material' / 1x)
C          Read elements in turn.  Echo and count.
      nelmts = 0
      do 25 i = 1,3
   25   nold(i) = 0
   30 read (*,*, end=911) nnew, srcnew, matnew
      if (nnew(1) .ne. nold(1) .or.
     *    nnew(2) .ne. nold(2) .or.
     *    nnew(3) .ne. nold(3)) then
	nelmts = nelmts + 1
	do 35 i = 1,3
	  nvtx(i,nelmts) = nnew(i)
   35     nold(i) = nnew(i)
	source(nelmts) = srcnew
	materl(nelmts) = matnew
	write (*, 1180) nnew, srcnew, matnew
	go to 30
      endif
 1180 format (1x, i3, 2i6, 2x, f10.5, i7)
C
C     Read list of fixed potential values and print.
C          Print header to start.
  120 write (*, 1200)
 1200 format (1x // 4x, 'Input fixed potentials' // 6x,
     1        'node', 12x, 'value' / 1x)
C          Declare all nodes to start off unconstrained.
      do 40 m = 1,nodes
	potent(m) = 0.
	constr(m) = .false.
   40   continue
C          Read and echo input.
      nconst = 0
      iold = 0
   60 read (*,*, end=911) inew, potnew
      if (inew .ne. iold) then
	nconst = nconst + 1
	constr(inew) = .true.
	potent(inew) = potnew
	iold = inew
	write (*, 1210) inew, potnew
	go to 60
      endif
 1210 format (6x, i3, 9x, f10.5)
C
C     Return to calling program.
  900 return
  911 call errexc('MESHIN', 1)
      end
C
C*************************************************************
C
      Subroutine output
C
C*************************************************************
C
C     Prints the results on the standard output device.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75)
      parameter (TINY = 1.e-35, HUGE = 1.e+35)
      logical constr
      common /problm/ nodes, nelmts,
     *  x(MAXNOD), y(MAXNOD), constr(MAXNOD), potent(MAXNOD),
     *  nvtx(3,MAXELM), source(MAXELM), materl(MAXELM)
      common /matrix/ p(MAXNOD,MAXNOD), resid(MAXNOD),
     *                step(MAXNOD), b2max
      common /workng/ pel(3,3), sel(3,3), tel(3,3), e(3)
C=============================================================
C
C     Print the nodes and the output potential values.
      write (*,1000) (i, x(i), y(i), potent(i),
     1                                 i = 1, nodes)
 1000 format (1x //// 12x, 'Final solution' // 3x, 'i', 8x,
     1        'x', 9x, 'y', 7x, 'potential' // (1x, i3, 2x,
     2        f10.5, F10.5, 3X, f10.5))
C
      return
      end
