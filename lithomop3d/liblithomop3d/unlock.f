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
      subroutine unlock(b,btot,id,iwink,idhist,ibond,bond,histry,
     & nstep,numnp,nwink,nhist,neq,numdif,lastep,iflag)
c
c       program to remove the winkler forces and leave the applied
c       forces at a specified step.  code uses load history factors
c       for forces, if they are defind
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
      integer nstep,numnp,nwink,nhist,neq,numdif,lastep,iflag
      integer id(ndof,numnp),iwink(2,nwink),idhist(numnp)
      integer ibond(ndof,numnp)
      double precision b(neq),btot(neq),bond(ndof,numnp)
      double precision histry(nhist,lastep+1)
c
c...  local constants
c
      double precision cutoff
      parameter(cutoff=0.8d0)
c
c...  local variables
c
      integer i,j,k,idk,mode,ihist,idforc
      double precision diff,forc
c
cdebug      write(6,*) "Hello from unlock_f!"
c
      do i=1,ndof
        do j=1,numnp
          if(id(i,j).ne.izero) then
            do k=1,nwink
              idk=iwink(2,k)
              if(id(i,j).eq.idk) then
                mode=iwink(1,k)
                if(mode.lt.izero) then
                  ihist=-mode
                  diff=histry(ihist,nstep+1)-histry(ihist,nstep)
                  if(diff.lt.-cutoff) then
                    if(iflag.eq.ione) then
                      idforc=ibond(i,j)
                      forc=bond(i,j)
                    else if(iflag.eq.itwo) then
                      idforc=idhist(j)
                      forc=zero
                      if(numdif.ne.izero) forc=bond(i,j)
                    end if
                    if(idforc.eq.izero) then
                      b(idk)=-btot(idk)+forc
                    else
                      b(idk)=-btot(idk)+forc*histry(idforc,nstep+1)
                    end if
                  end if
                end if
                go to 20
              end if
            end do
          end if
20        continue
        end do
      end do
      return
      end
c
c version
c $Id: unlock.f,v 1.2 2004/07/12 21:04:56 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
