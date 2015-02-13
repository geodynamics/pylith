// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file modulesrc/feassemble/IntegratorElasticity.i
 *
 * @brief Python interface to C++ abstract IntegratorElasticity object.
 */

namespace pylith {
  namespace feassemble {

    class IntegratorElasticity : public pylith::feassemble::Integrator
    { // IntegratorElasticity

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      IntegratorElasticity(void);

      /// Destructor
      virtual
      ~IntegratorElasticity(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set material.
       *
       * @param m Elastic material.
       */
      void material(pylith::materials::ElasticMaterial* m);
      
      /** Determine whether we need to recompute the Jacobian.
       *
       * @returns True if Jacobian needs to be recomputed, false otherwise.
       */
      bool needNewJacobian(void);
      
      /** Initialize integrator.
       *
       * @param mesh Finite-element mesh.
       */
      void initialize(const pylith::topology::Mesh& mesh);
      
      /** Update state variables as needed.
       *
       * @param t Current time
       * @param fields Solution fields
       * @param mesh Finite-element mesh
       */
      void updateStateVars(const PylithScalar t,
			   pylith::topology::SolutionFields* const fields);
      
      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      virtual
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of vertex field.
       * @param mesh Finite-element mesh for problem.
       * @param fields Fields manager.
       * @returns Cell field.
       */
      const pylith::topology::Field& cellField(const char* name,
					       const pylith::topology::Mesh& mesh,
					       pylith::topology::SolutionFields* const fields =0);
      
      /** Get output fields.
       *
       * @returns Output (buffer) fields.
       */
      const pylith::topology::Fields* outputFields(void) const;

    }; // IntegratorElasticity

  } // feassemble
} // pylith


// End of file 
