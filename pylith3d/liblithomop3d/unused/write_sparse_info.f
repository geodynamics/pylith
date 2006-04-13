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
      subroutine write_sparse_info(neq,nnz,nmin,nmax,wavg,idout,kw,
     & ofile)
c
c...subroutine to write global control parameters
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer neq,nnz,nmin,nmax,idout,kw
      double precision wavg
      character ofile*(*)
c
cdebug      write(6,*) "Hello from write_sparse_info_f!"
c
      if(idout.gt.0) then
        open(kw,file=ofile,status="old",access="append")
        write(kw,3000) neq,nnz-1,nmin,nmax,wavg
        close(kw)
      end if
c
      return
c
 3000 format(///,
     &' e q u a t i o n   s y s t e m   d a t a               ',   //,
     & 5x, 
     &' number of equations  . . . . . . . . . . . . . (neq) =',i8,//,
     & 5x,
     &' number of nonzero terms in stiffness matrix. . (nnz) =',i8,//,
     & 5x,
     &' minimum number of nonzero terms per row. . . .(nmin) =',i8,//,
     & 5x, 
     &' maximum number of nonzero terms per row. . . .(nmax) =',i8,//,
     & 5x, 
     &' average number of nonzero terms per row. . . .(wavg) =',f8.2)
c
      end
c
c version
c $Id: write_sparse_info.f,v 1.1 2005/08/05 19:58:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
