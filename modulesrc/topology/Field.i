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
 * @file modulesrc/topology/Field.hh
 *
 * @brief Python interface to C++ Field object.
 */

namespace pylith {
  namespace topology {

    class Field : public FieldBase
    { // Field

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param mesh Finite-element mesh.
       */
      Field(const pylith::topology::Mesh& mesh);

      /// Destructor.
      ~Field(void);
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);

      /** Get mesh associated with field.
       *
       * @returns Finite-element mesh.
       */
      const pylith::topology::Mesh& mesh(void) const;

      /** Set label for field.
       *
       * @param value Label for field.
       */
      void label(const char* value);

      /** Get label for field.
       *
       * @returns Label for field.
       */
      const char* label(void) const;
      
      /** Set vector field type
       *
       * @param value Type of vector field.
       */
      void vectorFieldType(const pylith::topology::FieldBase::VectorFieldEnum value);

      /** Get vector field type
       *
       * @returns Type of vector field.
       */
      pylith::topology::FieldBase::VectorFieldEnum vectorFieldType(void) const;

      /** Set scale for dimensionalizing field.
       *
       * @param value Scale associated with field.
       */
      void scale(const PylithScalar value);

      /** Get scale for dimensionalizing field.
       *
       * @returns Scale associated with field.
       */
      PylithScalar scale(void) const;
      
      /** Set flag indicating whether it is okay to dimensionalize field.
       *
       * @param value True if it is okay to dimensionalize field.
       */
      void dimensionalizeOkay(const bool value);
      
      /** Set flag indicating whether it is okay to dimensionalize field.
       *
       * @param value True if it is okay to dimensionalize field.
       */
      bool dimensionalizeOkay(void) const;
      
      /** Get spatial dimension of domain.
       *
       * @returns Spatial dimension of domain.
       */
      int spaceDim(void) const;
      
      /** Get the number of points in the chart.
       *
       * @returns the chart size.
       */
      int chartSize(void) const;
      
      /** Get the number of degrees of freedom.
       *
       * @returns the number of degrees of freedom.
       */
      int sectionSize(void) const;
      
      /** Has section been setup?
       *
       * @returns True if section has been setup.
       */
      bool hasSection(void) const;

      /// Set chart for solution.
      void setupSolnChart(void);
      
      /** Set default DOF for solution.
       *
       * @param fiberDim Total number of components in solution.
       * @param subfieldName Name of subfield for DOF.
       */
      void setupSolnDof(const int fiberDim,
			const char* subfieldName ="displacement");

      /** Create PETSc section and set chart and fiber dimesion.
       *
       * @param domain Type of points over which to define section.
       * @param dim Fiber dimension for section.
       * @param stratum Stratum depth (for vertices) and height (for cells).
       */
      void newSection(const pylith::topology::FieldBase::DomainEnum domain,
		      const int fiberDim,
		      const int stratum =0);

      /** Create section with same layout (fiber dimension and
       * constraints) as another section. This allows the layout data
       * structures to be reused across multiple fields, reducing memory
       * usage.
       *
       * @param sec Section defining layout.
       */
      void cloneSection(const Field& src);

      /** Add subfield to current field.
       *
       * Should be followed by calls to subfieldsSetup() and subfieldSetDof().
       *
       * @param name Name of subfield.
       * @param numComponents Number of components in subfield.
       * @param fieldType Type of vector field.
       * @param scale Scale for dimensionalizing field.
       */
      void subfieldAdd(const char *name, 
		       int numComponents,
		       const VectorFieldEnum fieldType,
		       const PylithScalar scale =1.0);

      /** Setup sections for subfields.
       *
       * Should be preceded by calls to subfieldAdd() and followed by calls to subfieldSetDof().
       */
      void subfieldsSetup(void);

      /** Convenience method for setting number of DOF (fiberdim) for subfield at points.
       *
       * Should be preceded by calls to subfieldAdd() and subfieldsSetup().
       *
       * @param name Name of subfield.
       * @param domain Point classification for subfield.
       * @param fiberDim Number of subfield components per point.
       */
      void subfieldSetDof(const char *name, 
			  const pylith::topology::FieldBase::DomainEnum domain, 
			  const int fiberDim);

      /// Clear variables associated with section.
      void clear(void);

      /// Allocate field.
      void allocate(void);
      
      /// Zero section values (excluding constrained DOF).
      void zero(void);
      
      /// Zero section values (including constrained DOF).
      void zeroAll(void);
      
      /// Complete section by assembling across processors.
      void complete(void);

      /** Copy field values and metadata.
       *
       * @param field Field to copy.
       */
      void copy(const Field& field);
      
      /** Copy subfield values and its metadata to field;
       *
       * @param field Field to copy from.
       * @param name Name of subfield to copy.
       */
      void copySubfield(const Field& field,
			const char* name);

      /** Add two fields, storing the result in one of the fields.
       *
       * @param field Field to add.
       */
      void add(const Field& field);

      /** Dimensionalize field. Throws runtime_error if field is not
       * allowed to be dimensionalized.
       */
      void dimensionalize(void);
      
      /** Print field to standard out.
       *
       * @param label Label for output.
       */
      void view(const char* label);

      /** Create PETSc vector scatter for field. This is used to transfer
       * information from the "global" PETSc vector view to the "local"
       * PETSc section view.
       *
       * @param mesh Mesh associated with scatter.
       * @param context Label for context associated with vector.
       */
      void createScatter(const pylith::topology::Mesh& mesh,
			 const char* context ="");

      /** Get PETSc vector associated with field.
       *
       * @param context Label for context associated with vector.
       * @returns PETSc vector.
       */
      PetscVec vector(const char* context ="");

      /** Get PETSc vector associated with field.
       *
       * @param context Label for context associated with vector.
       * @returns PETSc vector.
       */
      const PetscVec vector(const char* context ="") const;

      /** Scatter section information across processors to update the
       * global view of the field.
       *
       * @param context Label for context associated with vector.
       */
      void scatterLocalToGlobal(const char* context ="") const;
      
      /** Scatter section information across processors to update the
       * global view of the field.
       *
       * @param vector PETSc vector to update.
       * @param context Label for context associated with vector.
       */
      void scatterLocalToGlobal(const PetscVec vector,
				const char* context ="") const;
      
      /** Scatter global information across processors to update the local
       * view of the field.
       *
       * @param context Label for context associated with vector.
       */
      void scatterGlobalToLocal(const char* context ="") const;
      
      /** Scatter global information across processors to update the local
       * view of the field.
       *
       * @param vector PETSc vector used in update.
       * @param context Label for context associated with vector.
       */
      void scatterGlobalToLocal(const PetscVec vector,
				const char* context ="") const;

    }; // Field

  } // topology
} // pylith


// End of file 
