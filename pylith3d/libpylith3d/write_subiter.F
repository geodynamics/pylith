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
      subroutine write_subiter(nprevdflag,kw,idout,ofile)
c
c...subroutine to write parameters controlling the subiteration
c   process
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nprevdflag,kw,idout
      character ofile*(*)
c
c...  local variables
c
      if(idout.gt.0) then
        call open_append_stop(kw,ofile,"old","","write_subiter")
        write(kw,700) nprevdflag
        close(kw)
      end if
c
      return
c
700   format(///' l i n e a r   s o l u t i o n   ',
     & 'i n f o r m a t i o n'//,5x,
     & ' initial solution flag  . . . . . . . .(nprevdflag) =',i5,//,5x,
     & '    eq.0,  use zero initial guess                    ',/,5x,
     & '    eq.1,  use displacements from previous step      ')
      end
