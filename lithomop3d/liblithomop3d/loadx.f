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
      subroutine loadx(b,diforc,histry,idx,idhist,ndof,neq,numnp,
     & nhist,nstep,lastep,idout,kto,kw)
c
c...program to transfer differential forces into global load vector
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ndof,neq,numnp,nhist,nstep,lastep,idout,kto,kw
      integer idx(ndof,numnp),idhist(numnp)
      double precision b(neq),diforc(ndof,numnp),histry(nhist,lastep+1)
c
c...  local variables
c
      integer j,i,k,ihist
      double precision diff
c
cdebug      write(6,*) "Hello from loadx_f!"
c
      do j=1,numnp
        ihist=idhist(j)
        if((ihist.gt.nhist).or.(ihist.lt.0)) then
          if(idout.gt.1) write(kw,1000) ihist,j
          write(kto,1000) ihist,j
          stop
        end if
        do i=1,ndof
          k=idx(i,j)
          if(k.ne.0) then
            if((ihist.eq.0).and.(nstep.eq.0)) b(k)=b(k)+diforc(i,j)
            if(ihist.gt.0) then
              diff=histry(ihist,nstep+1)
              if(nstep.gt.0) diff=diff-histry(ihist,nstep)
              b(k)=b(k)+diforc(i,j)*diff
            end if
          end if
        end do
      end do
      return
1000  format(//' fatal slippery node error!'//
     & ' attempt to use undefined load history # ',i5,
     & ' for differential forces at node # ',i7)
      end
c
c version
c $Id: loadx.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
