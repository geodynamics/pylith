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
      subroutine assign_wink(winkdef,wink,iwinkdef,iwinkid,iwink,id,
     & numnp,nwink,nwinke)
c
c....program for reading and printing data on winkler restoring forces
c
c          winkdef(ndof,numnp) = values of winkler restoring spring
c                             constant, force(i,j)=-wink(i,j)*disp(i,j)
c
c          iwinkdef(ndof,numnp) = application mode:
c                               iwink = 0, no winkler forces,
c                               iwink = 1, applied throuthout computation
c                               iwink = -n, uses load history factor n
c
c          After assigning equation numbers, the following arrays are
c          returned:
c
c          wink(nwink)   = value of spring constant only for degrees of
c                          freedom with restoring forces.
c          wink(2,nwink) = application mode (1) and equation number (2)
c                          for each winkler restoring force.
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
      integer numnp,nwink,nwinke
      integer iwinkdef(ndof,nwinke),iwinkid(nwinke),iwink(2,nwink)
      integer id(ndof,numnp)
      double precision winkdef(ndof,nwinke),wink(nwink)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n,nwtot,nnz
c
c...  included variable definitions
c
c
c...  loop over Winkler entries
c
      nwtot=izero
      do i=1,nwinke
        nnz=izero
        n=iwinkid(i)
        do j=1,ndof
          if(iwinkdef(j,i).ne.izero) then
            nnz=nnz+1
            nwtot=nwtot+1
            iwink(1,nwtot)=iwinkdef(j,i)
            iwink(2,nwtot)=id(j,n)
            wink(nwtot)=winkdef(j,i)
          end if
        end do
      end do
c
      return
      end
c
c version
c $Id: assign_wink.f,v 1.1 2005/04/16 00:35:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
