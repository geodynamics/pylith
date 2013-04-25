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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/TractPerturbation.i
 *
 * @brief Python interface to C++ TractPerturbation object.
 */

namespace pylith {
  namespace faults {

    class TractPerturbation : public pylith::bc::TimeDependent
    { // class TractPerturbation

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      TractPerturbation(void);
      
      /// Destructor.
      virtual
      ~TractPerturbation(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
      
      /** Set label for traction perturbation.
       *
       * @param value Label.
       */
      void label(const char* value);
      
      /** Get parameter fields.
       *
       * @returns Parameter fields.
       */
      const pylith::topology::Fields<topology::Field<pylith::topology::SubMesh> >* parameterFields(void) const;
      
      /** Initialize slip time function.
       *
       * @param faultMesh Finite-element mesh of fault.
       * @param faultOrientation Orientation of fault.
       * @param normalizer Nondimensionalization of scales.
       */
      void initialize(const pylith::topology::SubMesh& faultMesh,
		      const pylith::topology::Field<pylith::topology::SubMesh>& faultOrientation,
		      const spatialdata::units::Nondimensional& normalizer);
      
      /** Calculate spatial and temporal variation of value.
       *
       * @param t Current time.
       */
      void calculate(const PylithScalar t);
      
      /** Determine if perturbation has a given parameter.
       *
       * @param name Name of parameter field.
       * @returns True if perturbation has parameter field, false otherwise.
       */
      bool hasParameter(const char* name) const;

      /** Get vertex field with traction perturbation information.
       *
       * @param name Name of field.
       * @param fields Solution fields.
       *
       * @returns Traction vector field.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>&
      vertexField(const char* name,
		  pylith::topology::SolutionFields* const fields =0);
      
      // PROTECTED METHODS //////////////////////////////////////////////
    protected :

      /** Get label of boundary condition surface.
       *
       * @returns Label of surface (from mesh generator).
       */
      const char* _getLabel(void) const;

    }; // class TractPerturbation
    
  } // faults
} // pylith


// End of file 
