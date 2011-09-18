// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

%include "petscconf.h"

#if defined(PETSC_USE_REAL_SINGLE)
#define PYLITH_USE_SCALAR_SINGLE
typedef float PylithScalar;
#else
#define PYLITH_USE_SCALAR_DOUBLE
typedef double PylithScalar;
#endif


// End of file
