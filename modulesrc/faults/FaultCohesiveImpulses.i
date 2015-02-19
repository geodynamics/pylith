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

/** @file modulesrc/faults/FaultCohesiveImpulses.i
 *
 * @brief Python interface to C++ FaultCohesiveImpulses object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveImpulses : public FaultCohesiveLagrange
    { // class FaultCohesiveImpulses

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveImpulses(void);

      /// Destructor.
      ~FaultCohesiveImpulses(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Sets the spatial database for amplitudes of the impulses.
       *
       * @param db spatial database for amplitudes of impulses.
       */
      void dbImpulseAmp(spatialdata::spatialdb::SpatialDB* db);
      
      /** Set indices of fault degrees of freedom associated with
       * impulses.
       *
       * @param flags Array of indices for degrees of freedom.
       * @param size Size of array
       */
      %apply(int* INPLACE_ARRAY1, int DIM1) {
	(const int* flags, 
	 const int size)
	  };
      void impulseDOF(const int* flags,
		      const int size);  
      %clear(const int* flags, const int size);
      
      /** Set threshold for nonzero impulse amplitude.
       *
       * @param value Threshold for detecting nonzero amplitude.
       */
      void threshold(const PylithScalar value);
      
      /** Get number of impulses.
       *
       * Multiply by number of components to get total number of impulses.
       *
       * @returns Number of points with impulses.
       */
      int numImpulses(void) const;
      
      /** Get number of components for impulses at each point.
       *
       * Multiply by number of components to get total number of impulses.
       *
       * @returns Number of points with impulses.
       */
      int numComponents(void) const;
      
      /** Initialize fault. Determine orientation and setup boundary
       * condition parameters.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Direction perpendicular to along-strike direction that is 
       *   not collinear with fault normal (usually "up" direction but could 
       *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
       */
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]);
      
      /** Integrate contributions to residual term (r) for operator that
       * do not require assembly across cells, vertices, or processors.
       *
       * @param residual Field containing values for residual
       * @param t Current time
       * @param fields Solution fields
       */
      void integrateResidual(const pylith::topology::Field& residual,
			     const PylithScalar t,
			     pylith::topology::SolutionFields* const fields);
      
      /** Get vertex field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Vertex field.
       */
      const pylith::topology::Field& vertexField(const char* name,
						 const pylith::topology::SolutionFields* fields =0);
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Cell field.
       */
      const pylith::topology::Field& cellField(const char* name,
					       const pylith::topology::SolutionFields* fields =0);

    }; // class FaultCohesiveImpulses

  } // faults
} // pylith


// End of file 
