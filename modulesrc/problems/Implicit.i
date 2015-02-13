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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/problems/Implicit.i
 *
 * @brief Python interface to C++ Implicit.
 */

namespace pylith {
  namespace problems {

    class Implicit : public Formulation
    { // Implicit

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
    public :

      /// Constructor
      Implicit(void);
      
      /// Destructor
      ~Implicit(void);

      /// Compute rate fields (velocity and/or acceleration) at time t.
      void calcRateFields(void);

    }; // Implicit

  } // problems
} // pylith


// End of file 
