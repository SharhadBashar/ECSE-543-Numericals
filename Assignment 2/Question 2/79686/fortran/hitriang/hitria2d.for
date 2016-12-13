C*************************************************************
C************                                     ************
C************  High-order finite element program  ************
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
C        nve     =  number of nodes in each element
C
C      Global matrices and right-hand side
C      common /matrix/
C        s       =  global S-matrix for whole problem
C        t       =  global T-matrix for whole problem
C        rthdsd  =  right-hand side of system of equations
C
C      Temporary working arrays and variables
C      common /workng/
C        sel     =  S for one element (working array)
C        tel     =  T for one element (working array)
C        intg    =  integer working array
C        intg0   =  integer working array
C
C      Predefined problem size parameters (array dimensions):
C        MAXNOD  =  maximum number of nodes in problem
C        MAXELM  =  maximum number of elements in problem
C        MAXNVE  =  maximum number of vertices in any element
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, MAXNVE = 15)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(MAXNOD),
     *        constr(MAXNOD), potent(MAXNOD),
     *        nvtx(MAXNVE,MAXELM), source(MAXELM), nve(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), t(MAXNOD,MAXNOD),
     *        rthdsd(MAXNOD)
      common /workng/ sel(MAXNVE,MAXNVE), tel(MAXNVE,MAXNVE),
     *        intg(MAXNVE), intg0(MAXNVE)
C=============================================================
C
C     Set up triangular finite element matrices.
      call gettri
C
C     Fetch input data from input file.
      call meshin
C
C        Set global matrices and right side to all zeros.
      call matini(s, MAXNOD, MAXNOD)
      call matini(t, MAXNOD, MAXNOD)
      call vecini(rthdsd, MAXNOD)
C
C     Assemble element matrix.  Proceed element by element.
      do 40 ie = 1,nelmts
C
C        Construct element s and t matrices
	if (nve(ie) .eq.  3  .or.  nve(ie) .eq.  6  .or.
     *      nve(ie) .eq. 10  .or.  nve(ie) .eq. 15)
     *      call trielm(ie)
C
C        Embed matrices in global s; augment right side:
	call elembd(ie)
   40   continue
C
C     Solve the assembled finite element equations
      call eqsolv (s, potent, rthdsd, nodes, MAXNOD)
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
C     Reads input data file in  three parts:  nodes, elements,
C     fixed potentials.   Each part is followed by a line con-
C     taining nothing but a slant character /  which serves as
C     a data terminator.   Data lines are in free format, with
C     comma or blank  as acceptable  separators between items.
C     Character  data  (element  types!)  must be  encased  in
C     apostrophes.
C
C     Nodes:                                    x, y   values.
C     Elements:       type, total nodes, source, node numbers.
C     Potentials:                    node number, fixed value.
C
C     All data are echoed as read but little checking is done
C     for validity.
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, MAXNVE = 15)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(MAXNOD),
     *        constr(MAXNOD), potent(MAXNOD),
     *        nvtx(MAXNVE,MAXELM), source(MAXELM), nve(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), t(MAXNOD,MAXNOD),
     *        rthdsd(MAXNOD)
      common /workng/ sel(MAXNVE,MAXNVE), tel(MAXNVE,MAXNVE),
     *        intg(MAXNVE), intg0(MAXNVE)
C=============================================================
C
      Parameter (HUGE = 1.e+35)
      logical done
C
C     Read in the node list and echo input lines.
C             Start by printing the heading for nodes.
      write (*, 200)
C
C          Read and echo until null list (line with / only)
      nodes = 0
      xi0 = HUGE
      yi0 = HUGE
   20 read (*,*) xi, yi
      if (xi .ne. xi0  .or.  yi .ne. yi0) then
	nodes = nodes + 1
	x(nodes) = xi
	y(nodes) = yi
	xi0 = xi
	yi0 = yi
	write (*, 210) nodes, xi, yi
	go to 20
      endif
C
C     Read in the element list and echo all input as received.
C          Print element heading to start.
      write (*, 220)
C
C          Read elements in turn.
      nelmts = 0
      do 30 i = 1,MAXNVE
   30   intg0(i) = 0
   40 read (*,*) nods, sourci,
     *                      (intg(j), j = 1,min0(nods,MAXNVE))
      done = .true.
      do 50 i = 1,nods
   50    done = done .and. (intg(i) .eq. intg0(i))
      if (.not. done) then
	nelmts = nelmts + 1
	source(nelmts) = sourci
	nve(nelmts) = nods
	write (*,230) nelmts, nods, sourci,
     *                                   (intg(i), i = 1,nods)
	do 60 i = 1,nods
	  nvtx(i,nelmts) = intg(i)
	  intg0(i) = intg(i)
   60     continue
	go to 40
      endif
C
C     Read list of fixed potential values and print.
C          Print header first, then declare all nodes free.
      write (*, 240)
      call vecini(potent, MAXNOD)
      do 70 m = 1,nodes
   70    constr(m) = .false.
C
C          Read and echo input.
      last = 0
   80 read (*,*, end=90) i, xi
      if (i .ne. last) then
	constr(i) = .true.
	potent(i) = xi
	last = i
	write (*, 250) i, xi
	go to 80
      endif
C
   90 return
C
  200 format (1x // 8x, 'Input node list' // 2x, 'Node', 6x, 'x',
     *        11x, 'y' / 1x)
  210 format (1x, i3, 2(2x, f10.5))
  220 format (1x // 9x, 'Input element list' //
     * '  No   Nodes  Source', '  Element node numbers' / 1x)
  230 format (1x, i3, 2x, i4, f10.5, 1x, 15i3)
  240 format (1x // 4x, 'Input fixed potentials' // 6x,
     *        'Node', 12x, 'Value' / 1x)
  250 format (6x, i3, 9x, f10.5)
C
      end
C
C*************************************************************
C
      Subroutine elembd(ie)
C
C*************************************************************
C
C     Embeds single-element  s and t matrices currently in sel
C     and tel (in common block /workng/) in the global  matrix
C     s.  Argument ie is the element number.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, MAXNVE = 15)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(MAXNOD),
     *        constr(MAXNOD), potent(MAXNOD),
     *        nvtx(MAXNVE,MAXELM), source(MAXELM), nve(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), t(MAXNOD,MAXNOD),
     *        rthdsd(MAXNOD)
      common /workng/ sel(MAXNVE,MAXNVE), tel(MAXNVE,MAXNVE),
     *        intg(MAXNVE), intg0(MAXNVE)
C=============================================================
C
C     Run through element  s and t matrices  sel and tel, aug-
C     menting s, t, and the right-hand side as appropriate.
C
      nvertc = nve(ie)
      do 60 i = 1,nvertc
	 irow = nvtx(i,ie)
C
C           Does row correspond to a fixed potential?
	 if (constr(irow)) go to 50
C
C           No, potential may vary.  Do all nvertc columns.
	 do 40 j = 1,nvertc
	    icol = nvtx(j,ie)
C
C           Does column correspond to a fixed potential?
	    if (constr(icol)) go to 30
C
C           No; so augment s, t and rthdsd.
	    s(irow,icol) = s(irow,icol) + sel(i,j)
	    t(irow,icol) = t(irow,icol) + tel(i,j)
	    rthdsd(irow) = rthdsd(irow) + tel(i,j) * source(ie)
	    go to 40
C
C           Yes; so augment right side only:
   30       continue
	    rthdsd(irow) = rthdsd(irow) + tel(i,j) * source(ie)
     1                            - sel(i,j) * potent(icol)
   40       continue
	 go to 60
C
C           Constrained row number.  Set global s and rthdsd.
   50    continue
	 s(irow,irow) = 1.
	 rthdsd(irow) = potent(irow)
C
   60    continue
C
C     All done -- return to calling program.
      return
      end
C
C*************************************************************
C
      Subroutine output
C
C*************************************************************
C
C     Outputs problem and results to standard output stream.
C
C=============================================================
C     Global declarations -- same in all program segments
C=============================================================
      parameter (MAXNOD = 50, MAXELM = 75, MAXNVE = 15)
      logical constr
      common /problm/ nodes, nelmts, x(MAXNOD), y(maxnod),
     *        constr(MAXNOD), potent(MAXNOD),
     *        nvtx(MAXNVE,MAXELM), source(MAXELM), nve(MAXELM)
      common /matrix/ s(MAXNOD,MAXNOD), t(MAXNOD,MAXNOD),
     *        rthdsd(MAXNOD)
      common /workng/ sel(MAXNVE,MAXNVE), tel(MAXNVE,MAXNVE),
     *        intg(MAXNVE), intg0(MAXNVE)
C=============================================================
C
C     Print the nodes and the output potential values.
C
      write (*, 1000) (i, x(i), y(i), potent(i),
     1                                 i = 1, nodes)
 1000 format (1x /// 12x, 'Final solution' // '   i', 8x, 'x',
     1        9x, 'y', 7x, 'Potential' // (1x, i3, 2x, f10.5,
     2        F10.5, 3x, f10.5))
C
      return
      end
