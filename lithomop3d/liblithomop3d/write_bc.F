c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c  All worldwide rights reserved.  A license to use, copy, modify and
c  distribute this software for non-commercial research purposes only
c  is hereby granted, provided that this copyright notice and
c  accompanying disclaimer is not modified or removed from the software.
c
c  DISCLAIMER:  The software is distributed "AS IS" without any express
c  or implied warranty, including but not limited to, any implied
c  warranties of merchantability or fitness for a particular purpose
c  or any warranty of non-infringement of any current or pending patent
c  rights.  The authors of the software make no representations about
c  the suitability of this software for any particular purpose.  The
c  entire risk as to the quality and performance of the software is with
c  the user.  Should the software prove defective, the user assumes the
c  cost of all necessary servicing, repair or correction.  In
c  particular, neither Rensselaer Polytechnic Institute, nor the authors
c  of the software are liable for any indirect, special, consequential,
c  or incidental damages related to the software, to the maximum extent
c  the law permits.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine write_bc(bond,ibond,numnp,kw,idout,ofile,ierr,errstrng)
c
c...  subroutine to write boundary conditions to ascii output file.
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file
c         4:  Write error
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numnp,kw,idout,ierr
      integer ibond(ndof,numnp)
      double precision bond(ndof,numnp)
      character ofile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  local constants
c
      character ltype(4)*4
      data ltype/'free','disp','velo','forc'/
c
c...  intrinsic functions
c
      intrinsic mod,index
c
c...  local variables
c
      integer ihist(3),itype(3)
      integer i,n,nlines,npage,imode
      logical nonzed
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
      if(idout.eq.izero) return
c
c...  output BC to ascii file, if requested
c
c...  at present, it is assumed that all f77 compilers support
c     access="append" and all f95 compilers support position="append".
c     If exceptions are found we will have to use FORTRAN_COMPILER_NAME
c     to differentiate.
c
#if defined FORTRAN_TYPE_F77
      open(kw,file=ofile,status="old",access="append",err=10)
#else
      open(kw,file=ofile,status="old",position="append",err=10)
#endif
      write(kw,1000) (labeld(i),i=1,ndof)
      write(kw,2000)
      nlines=0
      npage=50
      do n=1,numnp
        nonzed=.false.
        do i=1,ndof
          if(ibond(i,n).ne.0) nonzed=.true.
        end do
        if(nonzed) then
          nlines=nlines+1
          if(mod(nlines,npage).eq.0) then
            write(kw,1000) (labeld(i),i=1,ndof)
            write(kw,2000)
          end if
          do i=1,ndof
            imode=ibond(i,n)
            ihist(i)=imode/10
            itype(i)=imode-10*ihist(i)
          end do
          write(kw,3000,err=20) n,(bond(i,n),ltype(itype(i)+1),ihist(i),
     &     i=1,ndof)
        end if
      end do
      close(kw)
c
c...  normal return
c
      return
c
c...  error opening output file
c
10    continue
        ierr=2
        errstrng="write_bc"
        close(kw)
        return
c
c...  write error
c
20    continue
        ierr=4
        errstrng="write_bc"
        close(kw)
        return
c
1000  format(1x,///,'  n o d a l   f o r c e s   a n d   d i s p l a',
     1 ' c e m e n t s',//,
     2 '      key to boundary condition codes:',//,
     3 '          free = unconstrained degree of freedom; the',/,
     4 '                 numerical value is meaningless.',/,
     5 '          disp = fixed displacement',/,
     6 '          velo = constant velocity',/,
     7 '          forc = applied (constant) force',//,
     8 '      hfac = load history factor applied (0 if none)',///,
     9 ' node # ',3(4x,a4,6x,'type  hfac  '))
2000  format(' ')
3000  format(1x,i7,2x,3(1pe12.5,2x,a4,1x,i3,4x))
      end
c
c version
c $Id: write_bc.F,v 1.1 2005/08/05 19:58:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
