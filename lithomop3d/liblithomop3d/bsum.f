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
      subroutine bsum(bextern,btraction,bgravity,bconcforce,bintern,
     & bresid,nextflag,ntractflag,ngravflag,nconcflag,neq)
c
c...subroutine to sum all contributions to the external force vector,
c   then subtract the internal force vector from the external vector
c   to obtain the residual.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
      double precision alpha
      parameter(alpha=-1.0d0)
c
c...  subroutine arguments
c
      integer nextflag,ntractflag,ngravflag,nconcflag,neq
      double precision bextern(nextflag*neq),btraction(ntractflag*neq)
      double precision bgravity(ngravflag*neq),bconcforce(nconcflag*neq)
      double precision bintern(neq),bresid(neq)
c
cdebug      integer idb
c
cdebug      write(6,*) "Hello from bsum_f!"
c
      if(nextflag.ne.izero) then
        call fill(bextern,zero,neq)
        if(ntractflag.ne.izero) call daxpy(neq,one,btraction,ione,
     &   bextern,ione)
        if(ngravflag.ne.izero) call daxpy(neq,one,bgravity,ione,
     &   bextern,ione)
        if(nconcflag.ne.izero) call daxpy(neq,one,bconcforce,ione,
     &   bextern,ione)
        call dcopy(neq,bextern,ione,bresid,ione)
cdebug        write(6,*) "bextern:"
cdebug        write(6,*) (bextern(idb),idb=1,200)
cdebug        if(ntractflag.ne.izero) then
cdebug          write(6,*) "btraction:"
cdebug          write(6,*) (btraction(idb),idb=1,200)
cdebug        end if
cdebug        if(ngravflag.ne.izero) then
cdebug          write(6,*) "bgravity:"
cdebug          write(6,*) (bgravity(idb),idb=1,200)
cdebug        end if
cdebug        if(nconcflag.ne.izero) then
cdebug          write(6,*) "bconcforce:"
cdebug          write(6,*) (bconcforce(idb),idb=1,200)
cdebug        end if
      else
        call fill(bresid,zero,neq)
      end if
      call daxpy(neq,alpha,bintern,ione,bresid,ione)
cdebug      write(6,*) "bintern:"
cdebug      write(6,*) (bintern(idb),idb=1,200)
cdebug      write(6,*) "bresid:"
cdebug      write(6,*) (bresid(idb),idb=1,200)
      return
      end
c
c version
c $Id: bsum.f,v 1.2 2005/02/23 23:50:57 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
