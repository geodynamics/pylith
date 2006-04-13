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
c...  labels_def.inc:  A header that defines the contents of the
c     labels array.  This is a temporary solution to avoid the
c     problems of passing character arrays from python to f77.
c     At present, a header like this is required for each geometry type.
c     This file contains a data statement and so should be included
c     after all type statements, etc.
c
c...  arrays containing parameters
c
      data labels/"Strsxx-tot ","Strsyy-tot ","Strszz-tot ",
     &            "Strsxy-tot ","Strsyz-tot ","Strsxz-tot ",
     &            "Tstnxx-tot ","Tstnyy-tot ","Tstnzz-tot ",
     &            "Tstnxy-tot ","Tstnyz-tot ","Tstnxz-tot ",
     &            "Vstnxx-tot ","Vstnyy-tot ","Vstnzz-tot ",
     &            "Vstnxy-tot ","Vstnyz-tot ","Vstnxz-tot ",
     &            "Pstnxx-tot ","Pstnyy-tot ","Pstnzz-tot ",
     &            "Pstnxy-tot ","Pstnyz-tot ","Pstnxz-tot ",
     &            "Strsxx-inc ","Strsyy-inc ","Strszz-inc ",
     &            "Strsxy-inc ","Strsyz-inc ","Strsxz-inc ",
     &            "Tstnxx-inc ","Tstnyy-inc ","Tstnzz-inc ",
     &            "Tstnxy-inc ","Tstnyz-inc ","Tstnxz-inc ",
     &            "Vstnxx-inc ","Vstnyy-inc ","Vstnzz-inc ",
     &            "Vstnxy-inc ","Vstnyz-inc ","Vstnxz-inc ",
     &            "Pstnxx-inc ","Pstnyy-inc ","Pstnzz-inc ",
     &            "Pstnxy-inc ","Pstnyz-inc ","Pstnxz-inc ",
     &            "Strsxx-rate","Strsyy-rate","Strszz-rate",
     &            "Strsxy-rate","Strsyz-rate","Strsxz-rate",
     &            "Tstnxx-rate","Tstnyy-rate","Tstnzz-rate",
     &            "Tstnxy-rate","Tstnyz-rate","Tstnxz-rate",
     &            "Vstnxx-rate","Vstnyy-rate","Vstnzz-rate",
     &            "Vstnxy-rate","Vstnyz-rate","Vstnxz-rate",
     &            "Pstnxx-rate","Pstnyy-rate","Pstnzz-rate",
     &            "Pstnxy-rate","Pstnyz-rate","Pstnxz-rate"/
c
c...  scalar parameters associated with arrays
c
c
c version
c $Id: labels_def.inc,v 1.3 2005/03/23 17:16:09 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
