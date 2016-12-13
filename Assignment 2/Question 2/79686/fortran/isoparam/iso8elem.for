C*************************************************************
C
      Subroutine iso8el(ie)
C
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C
C     Constructs  the element matrices s, t  for an isoparame-
C     tric element with 8 nodes.  ie = element number.
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
      dimension alf(8), duv(2,8), tj(2,2), tji(2,2), tjti(2,2)
      dimension wq(5), sq(5)
C      common /gauss5/ wq, sq
      common /geom8/ alf, duv, tj, tji, tjti, detj
C     Gaussian quadrature formulae for m = 5
      data wq/0.2369269, 0.4786287, 0.5688889,
     &                                   0.4786287, 0.2369269/
      data sq/-0.9061799, -0.5384693, 0.0000000,
     &                                   0.5384693, 0.9061799/
C
C     Make sure there are 8 nodes
      if (nve(ie) .ne. 8)  call errexc('ISO8EL', nve(ie))
C
C     Clear arrays, then start the work
      do 30 i = 1,8
	do 20 j = 1,8
	  sel(i,j) = 0.
	  tel(i,j) = 0.
   20     continue
   30   continue
C
C     Compute element matrix entries sel(i,j), tel(i,j)
C
C     Indexing:  i counts matrix rows
C                j counts matrix columns
C                m counts quadrature nodes in u
C                n counts quadrature nodes in v
C
      do 280 m = 1,5
	do 270 n = 1,5
	  u = sq(m)
	  v = sq(n)
	  call i8geom(ie, u, v)
C           Compute sel, tel contributions
	  do 230 i = 1,nve(ie)
	    do 220 j = 1,i
	      tel(i,j) = tel(i,j) +
     &                       wq(m)*wq(n)*alf(i)*alf(j)*detj
	      sum = 0.
	      do 130 is = 1,2
		do 120 js = 1,2
		  sum = sum + duv(is,i) *
     &                        tjti(is,js) * duv(js,j)
  120             continue
  130           continue
	      sel(i,j) = sel(i,j) + wq(m)*wq(n)*sum*detj
  220         continue
  230       continue
  270     continue
  280   continue
C
C     Make sel, tel symmetric!
      do 310 i = 1,8
	do 300 j = 1,i
	 tel(j,i) = tel(i,j)
	 sel(j,i) = sel(i,j)
  300    continue
  310   continue
C
  900 return
      end
C
C*************************************************************
C
      Subroutine i8geom(ie, u, v)
C
C*************************************************************
C      Copyright (c) 1995  P.P. Silvester and R.L. Ferrari
C*************************************************************
C
C     Determines  geometric  parameters and  derivatives of an
C     isoparametric element with 8 nodes. ie = element number.
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
      dimension alf(8), duv(2,8), tj(2,2), tji(2,2), tjti(2,2)
      common /geom8/ alf, duv, tj, tji, tjti, detj
C
      Do 50 k = 1,8
C
C              Interpolation functions for 8-noded element
	alf(1) = (-1.+u+v)*(1.+u)*(1.+v)/4.
	alf(2) = (1.+u)*(1.+v)*(1.-u)/2
	alf(3) = (1.+v)*(1.-u)*(v-(1.+u))/4.
	alf(4) = (1.+v)*(1.-u)*(1.-v)/2
	alf(5) = -(1.+u+v)*(1.-u)*(1.-v)/4.
	alf(6) = (1.+u)*(1.-u)*(1.-v)/2
	alf(7) = (1.+u)*(1.-v)*(u-(1.+v))/4.
	alf(8) = (1.+u)*(1.+v)*(1.-v)/2.
C
C            Function derivatives with respect to u
	duv(1,1) = (2.*u+v*(1.+2.*u+v))/4.
	duv(1,2) = -u*(1.+v)
	duv(1,3) = (1.+v)*(2.*u-v)/4.
	duv(1,4) = (-1.+v)*(1.+v)/2
	duv(1,5) = (2.*u+v*(1.-(2.*u+v)))/4.
	duv(1,6) = u*(-1.+v)
	duv(1,7) = (2.*u+v*(v-(1.+2.*u)))/4.
	duv(1,8) = (1.+v)*(1.-v)/2.
C
C            Function derivatives with respect to v
	duv(2,1) = (1.+u)*(u+2.*v)/4.
	duv(2,2) = (1.+u)*(1.-u)/2
	duv(2,3) = (2.*v+u*(u-(1.+2.*v)))/4.
	duv(2,4) = v*(-1.+u)
	duv(2,5) = (1.-u)*(u+2.*v)/4.
	duv(2,6) = (-1.+u)*(1.+u)/2
	duv(2,7) = (1.+u)*(-u+2.*v)/4.
	duv(2,8) = -v*(1.+u)
   50   continue
C
C           Compute the Jacobian
      tj(1,1) = 0.
      tj(1,2) = 0.
      tj(2,1) = 0.
      tj(2,2) = 0.
      do 70 k = 1,8
	tj(1,1) = tj(1,1) + x(nvtx(k,ie))*duv(1,k)
	tj(1,2) = tj(1,2) + y(nvtx(k,ie))*duv(1,k)
	tj(2,1) = tj(2,1) + x(nvtx(k,ie))*duv(2,k)
	tj(2,2) = tj(2,2) + y(nvtx(k,ie))*duv(2,k)
   70   continue
C
C           Jacobian determinant, inverse, transpose-inverse
      detj = tj(1,1)*tj(2,2) - tj(1,2)*tj(2,1)
      tji(1,1) =  tj(2,2) / detj
      tji(1,2) = -tj(1,2) / detj
      tji(2,1) = -tj(2,1) / detj
      tji(2,2) =  tj(1,1) / detj
      do 100 is = 1,2
	do 90 js = 1,2
	  tjti(is,js) = 0.
	  do 80 ks = 1,2
	    tjti(is,js) = tjti(is,js) +
     &                     tji(ks,is) * tji(ks,js)
   80       continue
   90     continue
  100   continue
C
      return
      end
