      program test_makemsr
c
c...  simple driver to test makemsr routine given the output from
c     lnklst.
c
      implicit none
      integer neq,nnz,iwork
      parameter(neq=4521,nnz=296640,iwork=452100)
      integer ja(nnz),indx(neq),link(iwork),nbrs(iwork)
      integer nmin,nmax,i
      double precision wavg
c
      nmin=0
      nmax=0
      wavg=0.0d0
      open(15,file="makemsr.info")
      read(15,*) (indx(i),i=1,neq)
      read(15,*) (link(i),i=1,iwork)
      read(15,*) (nbrs(i),i=1,iwork)
      call makemsr(ja,indx,link,nbrs,neq,nnz,iwork,nmin,nmax,wavg)
      write(6,*) nmin,nmax,wavg
      stop
      end
c
c
c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
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
      subroutine makemsr(ja,indx,link,nbrs,neq,nnz,iwork,nmin,nmax,
     & wavg)
c
c      program to transform linked list into modified sparse row format
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer neq,nnz,iwork,nmin,nmax
      integer ja(nnz),indx(neq),link(iwork),nbrs(iwork)
      double precision wavg
c
c...  intrinsic functions
c
      intrinsic min,max,dble
c
c...  local variables
c
      integer i,loc,ncol
cdebug      integer idb
c
      write(6,*) "Hello from makemsr_f!"
c
      write(6,*) "neq,nnz,iwork,nmin,nmax:",neq,nnz,iwork,nmin,nmax
      call ifill(ja,izero,nnz)
cdebug      if(neq.ne.izero) stop
      ja(1)=neq+2
      nmin=ione
      nmax=ione
      wavg=neq
      do i=1,neq
        loc=indx(i)
        ncol=izero
 20     continue
        ja(ja(i)+ncol)=nbrs(loc)
        ncol=ncol+ione
        loc=link(loc)
        if(loc.gt.izero) goto 20
        nmin=min(nmin,ncol)
        nmax=max(nmax,ncol)
        wavg=wavg+dble(ncol)
        ja(i+1)=ja(i)+ncol
c
c...    sort entries in each row
c
        call isort(ja(i+1)-ja(i),ja(ja(i)))
      end do
      write(6,*) "After end of loop!"
      call flush(6)
cdebug      if(neq.ne.izero) stop
      nmin=nmin+1
      nmax=nmax+1
      if(neq.ne.izero) wavg=wavg/dble(neq)
      write(6,*) "nmin,nmax,wavg,neq:",nmin,nmax,wavg,neq
      call flush(6)
cdebug      write(6,*) "ja:",(ja(idb),idb=1,nnz)
      call flush(6)
cdebug      if(neq.ne.izero) stop
      return
      end
c
c version
c $Id: test_makemsr.f,v 1.1 2004/08/12 14:41:19 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
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
      subroutine isort(n,ia)
c
c...  subroutine to sort an integer array.  Adapted from Numerical
c     Recipes.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer n
      integer ia(n)
c
c... local variables
c
      integer l,ir,iia,j,i
c
      if(n.le.1) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          iia=ia(l)
        else
          iia=ia(ir)
          ia(ir)=ia(1)
          ir=ir-1
          if(ir.eq.1)then
            ia(1)=iia
            return
          end if
        end if
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ia(j).lt.ia(j+1))j=j+1
          end if
          if(iia.lt.ia(j))then
            ia(i)=ia(j)
            i=j
            j=j+j
          else
            j=ir+1
          end if
        go to 20
        end if
        ia(i)=iia
      go to 10
      end
c
c version
c $Id: test_makemsr.f,v 1.1 2004/08/12 14:41:19 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
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
      subroutine ifill(ia,ival,m)
c
c.... program to fill an interger array with the given value
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer m
      integer ia(m),ival
c
c...  local variables
c
      integer i
c
      do i=1,m
        ia(i)=ival
      end do
      return
      end
c
c version
c $Id: test_makemsr.f,v 1.1 2004/08/12 14:41:19 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
