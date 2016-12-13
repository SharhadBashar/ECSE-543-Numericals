C*************************************************************
C
      Subroutine trielm(ie)
C
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C
C     Constructs  the element  matrices s, t  for a triangular
C     element of order 1, 2, 3, or 4.   ie = element number.
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
C     High-order triangular finite element common blocks --
C         shared only by subroutines TRIELM and GETTRI
C=============================================================
      common /qtblok/ t1(6), t2(21), t3(55), t4(120),
     1                q1(6), q2(21), q3(55), q4(120)
C=============================================================
C
C       Set up corner indices for triangle;
      if      (nve(ie) .eq.  3) then
	norder = 1
      else if (nve(ie) .eq.  6) then
	norder = 2
      else if (nve(ie) .eq. 10) then
	norder = 3
      else if (nve(ie) .eq. 15) then
	norder = 4
      else
	call errexc('TRIELM', 1000 + ie)
      endif
      i = nvtx(1,ie)
      j = nvtx(nve(ie)-norder,ie)
      k = nvtx(nve(ie),ie)
C
C       Find element area and cotangents of included angles
      area = abs((x(j) - x(i)) * (y(k) - y(i)) -
     1           (x(k) - x(i)) * (y(j) - y(i))) / 2.
      if (area .le. 0.) then
	call errexc('TRIELM', 2000 + ie)
      endif
      ctng1 = ((x(j) - x(i)) * (x(k) - x(i)) +
     1         (y(j) - y(i)) * (y(k) - y(i))) / (2. * Area)
      ctng2 = ((x(k) - x(j)) * (x(i) - x(j)) +
     1         (y(k) - y(j)) * (y(i) - y(j))) / (2. * Area)
      ctng3 = ((x(i) - x(k)) * (x(j) - x(k)) +
     1         (y(i) - y(k)) * (y(j) - y(k))) / (2. * Area)
C
C     Compute element s and t matrix entries
      do 30 l = 1,nve(ie)
	do 20 m = 1,l
	  go to (1, 2, 3, 4), norder
    1     tel(l,m) = t1(locate(l,m)) * area
	  sel(l,m) = q1(locate(l,m)) * ctng1
     *         + q1(locate(indxp(l,1,1),indxp(m,1,1))) * ctng2
     *         + q1(locate(indxp(l,1,2),indxp(m,1,2))) * ctng3
	  go to 10
    2     tel(l,m) = t2(locate(l,m)) * area
	  sel(l,m) = q2(locate(l,m)) * ctng1
     *         + q2(locate(indxp(l,2,1),indxp(m,2,1))) * ctng2
     *         + q2(locate(indxp(l,2,2),indxp(m,2,2))) * ctng3
	  go to 10
    3     tel(l,m) = t3(locate(l,m)) * area
	  sel(l,m) = q3(locate(l,m)) * ctng1
     *         + q3(locate(indxp(l,3,1),indxp(m,3,1))) * ctng2
     *         + q3(locate(indxp(l,3,2),indxp(m,3,2))) * ctng3
	  go to 10
    4     tel(l,m) = t4(locate(l,m)) * area
	  sel(l,m) = q4(locate(l,m)) * ctng1
     *         + q4(locate(indxp(l,4,1),indxp(m,4,1))) * ctng2
     *         + q4(locate(indxp(l,4,2),indxp(m,4,2))) * ctng3
   10     tel(m,l) = tel(l,m)
	  sel(m,l) = sel(l,m)
   20     continue
   30   continue
   50 return
      end
C
C ************************************************************
C
      Function indxp (i, n, k)
C
C ************************************************************
C
C     Returns the  permuted row or column  index corresponding
C     to  row or column  i  of a  triangular element matrix of
C     order n, with the rotation permutation applied k times.
C
      dimension ntwist(83), jump(6)
C
      data ntwist/  3,  1,  2,  6,  3,  5,  1,  2,  4, 10,  6,
     *  9,  3,  5,  8,  1,  2,  4,  7, 15, 10, 14,  6,  9, 13,
     *  3,  5,  8, 12,  1,  2,  4,  7, 11, 21, 15, 20, 10, 14,
     * 19,  6,  9, 13, 18,  3,  5,  8, 12, 17,  1,  2,  4,  7,
     * 11, 16, 28, 21, 27, 15, 20, 26, 10, 14, 19, 25,  6,  9,
     * 13, 18, 24,  3,  5,  8, 12, 17, 23,  1,  2,  4,  7, 11,
     * 16, 22/
      data jump/0,  3,  9, 19, 34, 55/
C
C     "ntwist" contains row and column mappings  for the first
C     four orders of triangular element. "jump" contains poin-
C     ters to starting points in "ntwist" (minus 1).
C
      km3 = mod(k, 3)
      indxp = i
      if (km3 .eq. 0) go to 90
      do 70 l = 1, km3
	indxp = ntwist(jump(n)+indxp)
   70   continue
   90 continue
      return
      end
C
C ************************************************************
C
      Subroutine gettri
C
C ************************************************************
C
C     Sets up the universal element matrices qe and te, single
C     precision,  by fetching data  from files.  q1 and t1 are
C     the first order matrices,  q2 and t2 second order,  etc.
C     Common data storage area /qtblok/ is shared with TRIELM.
C=============================================================
C     High-order triangular finite element common blocks --
C         shared only by subroutines TRIELM and GETTRI
C=============================================================
      common /qtblok/ t1(6), t2(21), t3(55), t4(120),
     1                q1(6), q2(21), q3(55), q4(120)
C=============================================================
      Double precision value
C
C          Attach and read Q of order 1.
      open (Unit=10, err=911, file='QTRIA1.DAT')
      read (10,*) i, n
      do 10 k = 1,n
	read (10,*) i, j, junk, value
	q1(locate(i,j)) = value
   10   continue
      close(10)
C
C          Attach and read T of order 1.
      open (Unit=10, err=911, file='TTRIA1.DAT')
      read (10,*) i, n
      do 11 k = 10,n
	read (10,*) i, j, junk, value
	t1(locate(i,j)) = value
   11   continue
      close(10)
C
C          Attach and read Q of order 2.
      open (Unit=10, err=911, file='QTRIA2.DAT')
      read (10,*) i, n
      do 20 k = 1,n
	read (10,*) i, j, junk, value
	q2(locate(i,j)) = value
   20   continue
      close(10)
C
C          Attach and read T of order 2
      open (Unit=10, err=911, file='TTRIA2.DAT')
      read (10,*) i, n
      do 21 k = 1,n
	read (10,*) i, j, junk, value
	t2(locate(i,j)) = value
   21   continue
      close(10)
C
C          Attach and read Q of order 3.
      open (Unit=10, err=911, file='QTRIA3.DAT')
      read (10,*) i, n
      do 30 k = 1,n
	read (10,*) i, j, junk, value
	q3(locate(i,j)) = value
   30   continue
      close(10)
C
C          Attach and read T of order 3.
      open (Unit=10, err=911, file='TTRIA3.DAT')
      read (10,*) i, n
      do 31 k = 1,n
	read (10,*) i, j, junk, value
	t3(locate(i,j)) = value
   31   continue
      close(10)
C
C          Attach and read Q of order 4.
      open (Unit=10, err=911, file='QTRIA4.DAT')
      read (10,*) i, n
      do 40 k = 1,n
	read (10,*) i, j, junk, value
	q4(locate(i,j)) = value
   40   continue
      close(10)
C
C          Attach and read T of order 4
      open (Unit=10, err=911, file='TTRIA4.DAT')
      read (10,*) i, n
      do 41 k = 1,n
	read (10,*) i, j, junk, value
	t4(locate(i,j)) = value
   41   continue
      close(10)
  911 continue
      end
