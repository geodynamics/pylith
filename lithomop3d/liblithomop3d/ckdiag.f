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
      subroutine ckdiag(alnz,nzero,neq,nnz,idout,kto,kw)
c
c...program to test for zero or negative diagonals of the
c   stiffness matrix
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nzero,neq,nnz,idout,kto,kw
      double precision alnz(nnz)
c
c...  local variables
c
      integer iz,in,i
c
cdebug      write(6,*) "Hello from ckdiag_f!"
c
      iz=0
      in=0
      nzero=0
      do i=1,neq
        if(alnz(i).eq.zero) iz=iz+1
        if(alnz(i).lt.zero) in=in+1
      end do
      if((iz.ne.0).or.(in.ne.0)) then
        nzero=1
        if(idout.gt.1) write(kw,1000) iz,in
        write(kto,1000) iz,in
      end if
 1000 format(///' ***fatal error in stiffness matrix!'/
     & ' ',i7,' zero diagonals found'/
     & ' ',i7,' negative diagonals found'/)
      return
      end
c
c version
c $Id: ckdiag.f,v 1.2 2004/06/18 15:22:21 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
