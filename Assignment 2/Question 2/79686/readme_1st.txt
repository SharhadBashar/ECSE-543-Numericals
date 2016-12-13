Finite Elements for Electrical Engineers 

P.P. Silvester and R.L. Ferrari
Cambridge University Press (3rd edition)

This directory structure contains all the programs and all the
tabulated data given in the Third Edition. The programs differ
from those of the Second Edition: some slightly, some extensively,
and some are brand new. All have been tested with several of the
common Fortran 77 compilers and have also been checked with FTNCHEK
to verify syntax and semantics. 

The directory structure is specifically

 readme.1st 

 fortran/lossylin/*.*
   ,,   /simple2d/*.*
   ,,   /hitriang/*.*
   ,,   /saturab/*.*
   ,,   /isoparam/*.*
   ,,   /utility/*.*

 execut/*.*
 simplex/*.*

The individual files are all uncompressed and, excepting those in execut, 
are plain text.  The ftp procedure can only retrieve the files in one
directory at a time, and to download all of the files here takes some  
time.  If you have PKUNZIP 2.0 available, it might be better to go to 
the alternative ftp site ftp.cup.cam.ac.uk, where in the directory 
 
 pub/science/FiniteElements 

you will find the compressed file 

 fe4ee3.zip

This binary file can be transferred in one operation and subsequently 
be unzipped with the command

  pkunzip -d fe4ee3.zip <path>

to reproduce the whole directory structure above in one operation.

However, the procedure for transferring the uncompressed files here 
is as follows:

Create a disk for your local pc with the empty directories above.  
In the club.eng ftp navigate to one of the directories, say by means 
of the ftp command 

 cd fortran/lossylin

Change to the corresponding local directory by means of 

 lcd A:\fortran\lossylin

Issue the ftp retrieval command 

 mget *.* 

This will fill up the current local directory with the 
appropriate files.  Navigate your way around the the rest of the    
directories to fill up the empty local ones appropriately. (note  
that the execut/*.* files are binary)

PPS & RLF 5/9/96


