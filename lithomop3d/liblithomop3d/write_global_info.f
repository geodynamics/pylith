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
      subroutine write_global_info(title,idout,idsk,numnp,
     & icode,idebug,kw,kp,ofile,pfile)
c
c...subroutine to write global control parameters
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer idout,idsk,numnp,icode,idebug,kw,kp
      character title*(*),ofile*(*),pfile*(*)
c
      if(idout.gt.izero) then
        open(kw,file=ofile,status="new",access="append")
        write(kw,500)
        write(kw,1000) title,idout,idsk,ngem,numnp,nsd,ndof,icode,idebug
        close(kw)
      end if
      if(idsk.eq.ione) then
        open(kp,file=pfile,status="new",access="append")
        write(kp,2000) title
        write(kp,3000) ngem,numnp,nsd,ndof,nstr
      else if(idsk.eq.itwo) then
        open(kp,file=pfile,status="new",access="append",
     &   form="unformatted")
        write(kp) title
        write(kp) ngem,numnp,nsd,ndof,nstr
      end if
      close(kp)
c
      return
c
 500  format(//,7x,
     x'******************************************************',/,7x,
     x'*                                                    *',/,7x,
     x'*              LITHOMOP 0.7 OUTPUT                   *',/,7x,
     x'*                                                    *',/,7x,
     x'*        Copyright 2004 by Charles A. Williams.      *',/,7x,
     x'*          Rensselaer Polytechnic Institute.         *',/,7x,
     x'*                All rights reserved.                *',/,7x,
     x'*                                                    *',/,7x,
     x'******************************************************',//)
1000  format(1x,a70,///,
     x' c o n t r o l   i n f o r m a t i o n                 ',  //,5x,
     x' ascii file output option . . . . . . . . . . (idout) =',i5,//,
     x 5x,
     x'    eq.0,  no ascii output file                        ',   /,5x,
     x'    eq.1,  echo input                                  ',   /,5x,
     x'    eq.2,  full output                                 ',  //,5x,
     x' fedsk.dat output option. . . . . . . . . . . (idsk ) =',i5,//,
     x 5x,
     x'    eq.0,  no output                                   ',   /,5x,
     x'    eq.1,  ascii output                                ',   /,5x,
     x'    eq.2,  binary output                               ',  //,5x,
     x' geometry type  . . . . . . . . . . . . . . . (ngem ) =',i5,//,
     x 5x,
     x'    eq.0,  axisymmetric                                ',  / ,5x,
     x'    eq.1,  plane strain                                ',  / ,5x,
     x'    eq.2,  plane stress                                ',  / ,5x,
     x'    eq.3,  out-of-plane                                ',  / ,5x,
     x'    eq.4,  three-dimensional                           ',  //,5x,
     x' number of nodal points . . . . . . . . . . . (numnp) =',i7,//,
     x 5x,
     x' number of space dimensions . . . . . . . . . (nsd  ) =',i5,//,
     x 5x,
     x' number of degrees of freedom per node  . . . (ndof ) =',i5,//,
     x 5x,
     x' analysis code  . . . . . . . . . . . . . . . (icode) =',i5,/,5x,
     x'    eq.0,  set up storage, echo input data             ',  / ,5x,
     x'    eq.1,  factor stiffness, print diagonals           ',  / ,5x,
     x'    eq.2,  elastic solution only                       ',  / ,5x,
     x'    eq.3,  full viscoelastic solution                  ',  //,5x,
     x' level of debugging output . . . . . . . . . (idebug) =',i5,/,5x,
     x'    eq.0,  no special debug output                     ',  / ,5x,
     x'    eq.1,  print local and global force vectors        ',  //,5x,
     x/)
2000  format(a70)
3000  format(16i7)
c
      end
c
c version
c $Id: write_global_info.f,v 1.4 2004/11/01 20:42:04 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
