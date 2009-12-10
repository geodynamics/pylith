// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/feassemble/IntegratorElasticityLgDeform.i
 *
 * @brief Python interface to C++ abstract IntegratorElasticityLgDeform object.
 */

namespace pylith {
  namespace feassemble {

    class IntegratorElasticityLgDeform : 
      public pylith::feassemble::IntegratorElasticity
    { // IntegratorElasticityLgDeform

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :
      
      /// Constructor
      IntegratorElasticityLgDeform(void);

      /// Destructor
      virtual
      ~IntegratorElasticityLgDeform(void);
      
      /** Determine whether we need to recompute the Jacobian.
       *
       * @returns True if Jacobian needs to be recomputed, false otherwise.
       */
      bool needNewJacobian(void);
      
      /** Update state variables as needed.
       *
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      void updateStateVars(const double t,
			   pylith::topology::SolutionFields* const fields);
      
    }; // IntegratorElasticityLgDeform

  } // feassemble
} // pylith


// End of file 
