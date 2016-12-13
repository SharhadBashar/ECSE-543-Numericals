C***********************************************************
C     This demonstration program reads lists of integers,
C     separated by the  Fortran 77 null-read character \.
C     List length must not exceed LONG but zero is OK.
C***********************************************************
      Parameter (NULL=0, LONG=4)
      Dimension list(LONG)
   10 continue
C
      i = 0
      iold = NULL
   20 read (*,*, end=90) inew
      if (inew .ne. iold) then
	write (*,*) inew
	i = i + 1
	list(i) = inew
	iold = inew
	go to 20
      else
	write (*,100) i
  100   format (' Length =', i2, '. Next list please!')
      endif
C
      go to 10
   90 stop
      end
