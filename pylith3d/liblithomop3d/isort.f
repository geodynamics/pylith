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
c $Id: isort.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
