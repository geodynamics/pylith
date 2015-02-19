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

/** @file modulesrc/faults/FaultCohesiveKin.i
 *
 * @brief Python interface to C++ FaultCohesiveKin object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesiveKin : public FaultCohesiveLagrange
    { // class FaultCohesiveKin

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesiveKin(void);
      
      /// Destructor.
      virtual
      ~FaultCohesiveKin(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set kinematic earthquake sources.
       *
       * @param names Array of kinematic earthquake source names.
       * @param numNames Number of earthquake sources.
       * @param sources Array of kinematic earthquake sources.
       * @param numSources Number of earthquake sources.
       */
      %apply(const char* const* string_list, const int list_len){
	(const char* const* names,
	 const int numNames)
	  };
      void eqsrcs(const char* const* names,
		  const int numNames,
		  EqKinSrc** sources,
		  const int numSources);
      %clear(const char* const* names, const int numNames);
      
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

    }; // class FaultCohesiveKin

  } // faults
} // pylith


// End of file 
