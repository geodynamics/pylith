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
      subroutine convert_case(string,iopt)
c
c...  routine to convert a character string from upper to lower case,
c     or vice-versa, depending on iopt.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iopt
      character *(*) string
c
c...  defined constants
c
      character upper*26,lower*26
      data lower/"abcdefghijklmnopqrstuvwxyz"/
      data upper/"ABCDEFGHIJKLMNOPQRSTUVWXYZ"/
c
c...  intrinsic functions
c
      intrinsic len
c
c...  local variables
c
      integer i,j,lenstr
c
      lenstr=len(string)
      if(iopt.eq.1) then
        do i=1,lenstr
          do j=1,26
            if(string(i:i).eq.upper(j:j)) then
              string(i:i)=lower(j:j)
              go to 10
            end if
          end do
 10       continue
        end do
      else
        do i=1,lenstr
          do j=1,26
            if(string(i:i).eq.lower(j:j)) then
              string(i:i)=upper(j:j)
              go to 20
            end if
          end do
 20       continue
        end do
      end if
c
      return
      end
c
c version
c $Id: convert_case.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
