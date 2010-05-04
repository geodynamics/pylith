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

      /// Compute rate fields (velocity and/or acceleration) at time t.
      void calcRateFields(void);

      // PROTECTED MEMBERS //////////////////////////////////////////////
      protected :

      /// Setup rate fields.
      void _setupRateFields(void);

    }; // Implicit

  } // problems
} // pylith


// End of file 
