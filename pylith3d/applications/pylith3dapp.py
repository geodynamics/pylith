#!@INTERPRETER@
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2004  All Rights Reserved
#
#  Copyright 2004 Rensselaer Polytechnic Institute.
#  All worldwide rights reserved.  A license to use, copy, modify and
#  distribute this software for non-commercial research purposes only
#  is hereby granted, provided that this copyright notice and
#  accompanying disclaimer is not modified or removed from the software.
#
#  DISCLAIMER:  The software is distributed "AS IS" without any express
#  or implied warranty, including but not limited to, any implied
#  warranties of merchantability or fitness for a particular purpose
#  or any warranty of non-infringement of any current or pending patent
#  rights.  The authors of the software make no representations about
#  the suitability of this software for any particular purpose.  The
#  entire risk as to the quality and performance of the software is with
#  the user.  Should the software prove defective, the user assumes the
#  cost of all necessary servicing, repair or correction.  In
#  particular, neither Rensselaer Polytechnic Institute, nor the authors
#  of the software are liable for any indirect, special, consequential,
#  or incidental damages related to the software, to the maximum extent
#  the law permits.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# 


# main

if __name__ == "__main__":
    
    # re-create the PYTHONPATH at 'configure' time
    import sys
    path = '@PYTHONPATH@'.split(':')
    path.reverse()
    for dir in path:
        if dir:
            sys.path.insert(1, dir)

    # if we are embedding, insert the extension module in the
    # 'pylith3d' package
    try:
        import builtin_pylith3d
        sys.modules['pylith3d.pylith3d'] = builtin_pylith3d
    except ImportError:
        pass
    
    from pylith3d.Application import Application

    app = Application()
    app.run()
    

# version
__id__ = "$Id: pylith3dapp.py,v 1.2 2005/03/11 02:30:46 willic3 Exp $"

#  End of file 
