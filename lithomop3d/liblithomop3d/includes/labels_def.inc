c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2003  All Rights Reserved
c
c <LicenseText>
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  labels_def.inc:  A header that defines the contents of the
c     labels array.  This is a temporary solution to avoid the
c     problems of passing character arrays from python to f77.
c     At present, a header like this is required for each geometry type.
c     This file contains a data statement and so should be included
c     after all type statements, etc.
c
c...  arrays containing parameters
c
      data labels/"Strsxx","Strsyy","Strszz","Strsxy","Strsyz","Strsxz",
     &            "Tstnxx","Tstnyy","Tstnzz","Tstnxy","Tstnyz","Tstnxz",
     &            "Vstnxx","Vstnyy","Vstnzz","Vstnxy","Vstnyz","Vstnxz",
     &            "Pstnxx","Pstnyy","Pstnzz","Pstnxy","Pstnyz","Pstnxz"/
c
c...  scalar parameters associated with arrays
c
c
c version
c $Id: labels_def.inc,v 1.2 2005/03/22 23:24:23 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
