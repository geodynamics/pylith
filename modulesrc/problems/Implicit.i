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

    // PROTECTED METHODS ////////////////////////////////////////////////
    protected :

      /// Compute velocity at time t.
      void _calcVelocity(void);

    }; // Implicit

  } // problems
} // pylith


// End of file 
