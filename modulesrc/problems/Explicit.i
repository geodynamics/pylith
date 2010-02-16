// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file modulesrc/problems/Explicit.i
 *
 * @brief Python interface to C++ Explicit.
 */

namespace pylith {
  namespace problems {

    class Explicit : public Formulation
    { // Explicit

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
    public :

      /// Constructor
      Explicit(void);
      
      /// Destructor
      ~Explicit(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
    protected :

      /// Compute velocity at time t.
      void _calcVelocity(void);

    }; // Explicit

  } // problems
} // pylith


// End of file 
