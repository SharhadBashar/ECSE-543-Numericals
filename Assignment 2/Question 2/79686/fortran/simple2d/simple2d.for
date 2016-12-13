C*************************************************************
C************                                     ************
C************  First-order demonstration program  ************
C************                                     ************
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C
C      The subroutines  that make up this program  communicate
C      via named common blocks.  The variables in commons are:
C
C      Problem definition and solution
C      common /problm/
C        nodes   =  number of nodes used in problem
C        nelmts  =  number of elements in model
C        x, y    =  nodal coordinates
C        constr  =  logical, .true. for fixed potentials
C        potent  =  nodal potential array
C        nvtx    =  list of nodes for each element
C        source  =  source density in each element
C
C      Global matrix and right-hand side
C      common /matrix/
C        s       =  global s-matrix for whole problem
C        rthdsd  =  right-hand side of system of equations
C
C      Temporary working arrays and variables
C      common /workng/
C        sel     =  S for one element (working array)
C        tel     =  T for one element (working array)
C        intg    =  integer working array
C
C      Predefined problem size parameters (array dimensions):
C        MAXNOD  =  maximum number of nodes in problem
C        MAXELM  =  maximum number of elements in problem
C        HUGE    =  upper bound for any numeric value
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, HUGE = 1.E+35)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(MAXNOD),
     *                constr(MAXNOD), potent(MAXNOD),
     *                nvtx(3,MAXELM), source(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), rthdsd(MAXNOD)
      common /workng/ sel(3,3), tel(3,3), intg(3)
C=============================================================
C
C     Fetch input data from input file
      call meshin
C
C     Set global s-matrix and right side to all zeros.
      call matini(s, MAXNOD, MAXNOD)
      call vecini(rthdsd, MAXNOD)
C
C     Assemble global matrix, element by element.
      do 40 i = 1,nelmts
C           Construct element s and t matrices
	ie = i
	call elmatr(ie)
C          Embed matrices in global s; augment right side:
	call elembd(ie)
   40   continue
C
C     Solve the assembled finite element equations
      call eqsolv(s, potent, rthdsd, nodes, MAXNOD)
C
C     Print out the resulting potential values
      call output
C
      stop
      end
C
C*************************************************************
C
      Subroutine meshin
C
C*************************************************************
C
C     Reads input data file in three parts:   nodes, elements,
C     fixed potentials.  Each part is concluded by a line that
C     contains only the / character at its leftmost position.
C
C     Nodes:              node number, x, y.
C     Elements:           node numbers, source density.
C     Potentials:         node number, fixed value.
C
C     All data are echoed as read  but little checking is done
C     for validity.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, HUGE = 1.E+35)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(MAXNOD),
     *                constr(MAXNOD), potent(MAXNOD),
     *                nvtx(3,MAXELM), source(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), rthdsd(MAXNOD)
      common /workng/ sel(3,3), tel(3,3), intg(3)
C=============================================================
C
      Dimension nold(3), nnew(3)
C
C     Read in the node list and echo input lines.
C        Start by printing a heading for the node list.
      write (*, 1140)
 1140 format (1x // 8x, 'Input node list' / 3x, 'n', 8x, 'x',
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
 1160 format (1x // 6x, 'Input element list' / 3x, 'i', 5x,
     1        'j', 5x, 'k', 6x, 'Source' / 1x)
C
C          Read elements in turn.  Echo and count.
      nelmts = 0
      do 25 i = 1,3
   25   nold(i) = 0
   30 read (*,*, end=911) nnew, srcnew
      if (nnew(1) .ne. nold(1) .or.
     *    nnew(2) .ne. nold(2) .or.
     *    nnew(3) .ne. nold(3)) then
	nelmts = nelmts + 1
	do 35 i = 1,3
	  nvtx(i,nelmts) = nnew(i)
   35     nold(i) = nnew(i)
	source(nelmts) = srcnew
	write (*, 1180) nnew, srcnew
	go to 30
      endif
 1180 format (1x, i3, 2i6, 2x, g10.5)
C
C     Read list of fixed potential values and print.
C          Print header to start.
  120 write (*, 1200)
 1200 format (1x // 5x, 'Input fixed potentials' / 6x,
     1        'node', 12x, 'value' / 1x)
C          Declare all nodes to start off unconstrained.
      do 40 m = 1,nodes
	constr(m) = .false.
   40   continue
      call vecini(potent, nodes)
C
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
      return
  911 call errexc('MESHIN', 3)
      end
C
C*************************************************************
C
      Subroutine elmatr(ie)
C
C*************************************************************
C
C     Constructs element matrices  S and T for a single first-
C     order triangular finite element.  ie = element number.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, HUGE = 1.E+35)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(MAXNOD),
     *                constr(MAXNOD), potent(MAXNOD),
     *                nvtx(3,MAXELM), source(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), rthdsd(MAXNOD)
      common /workng/ sel(3,3), tel(3,3), intg(3)
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
   10       tel(l,m) = area / 12.
   20    tel(l,l) = 2. * tel(l,l)
C
C     Compute element S-matrix
      i1 = 1
      i2 = 2
      i3 = 3
      call matini(sel, 3, 3)
      do 50 nvrtex = 1,3
	 ctng = ((x(j) - x(i)) * (x(k) - x(i)) +
     1        (y(j) - y(i)) * (y(k) - y(i))) / (2. * area)
	 ctng2 = ctng / 2.
C
	 sel(i2,i2) = sel(i2,i2) + ctng2
	 sel(i2,i3) = sel(i2,i3) - ctng2
	 sel(i3,i2) = sel(i3,i2) - ctng2
	 sel(i3,i3) = sel(i3,i3) + ctng2
C
C          Permute row and column indices once
	 i4 = i1
	 i1 = i2
	 i2 = i3
	 i3 = i4
	 l = i
	 i = j
	 j = k
	 k = l
   50    continue
C
      return
      end
C
C*************************************************************
C
      Subroutine elembd(ie)
C
C*************************************************************
C
C     Embeds single-element S and T matrices currently in sel
C     and tel (in common block "workng") in the global matrix
C     s.  Argument ie is the element number.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, HUGE = 1.E+35)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(MAXNOD),
     *                constr(MAXNOD), potent(MAXNOD),
     *                nvtx(3,MAXELM), source(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), rthdsd(MAXNOD)
      common /workng/ sel(3,3), tel(3,3), intg(3)
C=============================================================
C
C     Run through element S and T matrices (sel and tel),
C     augmenting the global S and the right-hand side.
      do 60 i = 1,3
	irow = nvtx(i,ie)
C
C         Does row correspond to a fixed potential?
	if (constr(irow)) then
C           Constrained row number.  Set global s and rthdsd.
	  s(irow,irow) = 1.
	  rthdsd(irow) = potent(irow)
	else
C           No, potential is free to vary.  Do all 3 columns.
	  do 40 j = 1,3
	    icol = nvtx(j,ie)
C
C             Does column correspond to a fixed potential?
	    if (constr(icol)) then
C               Yes; so augment right side only:
	      rthdsd(irow) = rthdsd(irow) + tel(i,j)
     1                  * source(ie) - sel(i,j) * potent(icol)
	    else
C               No; so augment s and rthdsd.
	      s(irow,icol) = s(irow,icol) + sel(i,j)
	      rthdsd(irow) = rthdsd(irow)
     1                       + tel(i,j) * source(ie)
	    endif
   40     continue
	endif
   60   continue
C
C     All done -- return to calling program.
      return
      end
C
C
C*************************************************************
C
      Subroutine output
C
C*************************************************************
C
C     Prints the results on the standard output stream.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, HUGE = 1.E+35)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(MAXNOD),
     *                constr(MAXNOD), potent(MAXNOD),
     *                nvtx(3,MAXELM), source(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), rthdsd(MAXNOD)
      common /workng/ sel(3,3), tel(3,3), intg(3)
C=============================================================
C
C     Print the nodes and the output potential values.
C
      write (*,1000) (i, x(i), y(i), potent(i),
     1                                 i = 1, nodes)
 1000 format (1x // 12x, 'Final solution' / 3x, 'i', 8x,
     1        'x', 9x, 'y', 7x, 'potential' // (1x, i3, 2x,
     2        f10.5, F10.5, 3X, f10.5))
C
      return
      end
C
