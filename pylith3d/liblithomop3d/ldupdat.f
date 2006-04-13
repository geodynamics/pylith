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
      subroutine ldupdat(d,dx,tfault,dl,xl,ien,lmx,lmf,
     & nen,numnp,numfn,numslp)
c
c...subroutine to update local coordinates for large deformation
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
      integer nen,numnp,numfn,numslp
      integer ien(nen),lmx(ndof,nen),lmf(nen)
      double precision d(ndof,numnp),dx(ndof,numnp),tfault(ndof,numfn)
      double precision dl(ndof*nen),xl(nsd*nen)
c
c...  local variables
c
c
cdebug      write(6,*) "Hello from ldupdat_f"
c
      call ldisp(dl,d,ien,nen,numnp)
      if(numfn.ne.0) call adfldp(dl,lmf,tfault,nen,numfn)
      if(numslp.ne.0) call addsn(dl,dx,ien,lmx,nen,numnp)
      call daxpy(nen*nsd,one,dl,ione,xl,ione)
clater      if(iopt.eq.1.and.ldtmp.ne.0) then
clater        if(nsd.eq.ndof) then
clater          call daxpy(nen*nsd,one,dl,ione,xl,ione)
c
c...  Ugly (temporary) kludge to handle case where ndof.ne.nsd.
c     I could use this same setup for all cases, but it's probably
c     slower, and the special case only applies to out-of-plane
c     problems.
c
clater        else
clater          do i=1,nsd
clater            call daxpy(nen,one,dl(i),ndof,xl(i),nsd)
clater          end do
clater        end if
clater      end if
      return
      end
c
c version
c $Id: ldupdat.f,v 1.3 2004/08/12 01:34:40 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
