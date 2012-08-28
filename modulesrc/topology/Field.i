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
// Copyright (c) 2010-2012 University of California, Davis
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

// Typemaps for returning reference in operator+=. The default
// behavior is that the Python object will gain ownership. We want the
// C++ object to retain ownership. This could be alleviated by using a
// shared pointer.
%typemap(out) pylith::topology::Field<pylith::topology::Mesh>& {
  $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), $descriptor(pylith::topology::Field<pylith::topology::Mesh>*), 0 | 0);
  return $result;
 }
%typemap(out) pylith::topology::Field<pylith::topology::SubMesh>& {
  $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), $descriptor(pylith::topology::Field<pylith::topology::SubMesh>*), 0 | 0);
  return $result;
 }


namespace pylith {
  namespace topology {

    template<typename mesh_type>
    class Field : public FieldBase
    { // Field

      // PRIVATE TYPEDEFS ///////////////////////////////////////////////
    private:

      // Convenience typedefs
      typedef typename mesh_type::RealSection RealSection;
      typedef typename mesh_type::SieveMesh SieveMesh;
      typedef typename SieveMesh::label_sequence label_sequence;
      typedef typename RealSection::chart_type chart_type;
      
      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param mesh Finite-element mesh.
       */
      Field(const mesh_type& mesh);

      /// Destructor.
      ~Field(void);
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);

      /** Get mesh associated with field.
       *
       * @returns Finite-element mesh.
       */
      const mesh_type& mesh(void) const;

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
      void addDimensionOkay(const bool value);
      
      /** Set flag indicating whether it is okay to dimensionalize field.
       *
       * @param value True if it is okay to dimensionalize field.
       */
      bool addDimensionOkay(void) const;
      
      /** Get spatial dimension of domain.
       *
       * @returns Spatial dimension of domain.
       */
      int spaceDim(void) const;
      
      /** Get the number of sieve points in the chart.
       *
       * @returns the chart size.
       */
      int chartSize(void) const;
      
      /** Get the number of degrees of freedom.
       *
       * @returns the number of degrees of freedom.
       */
      int sectionSize(void) const;
      
      /// Create sieve section.
      void newSection(void);

      /** Create sieve section and set chart and fiber dimesion.
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

      void addField(const char *name, int numComponents);

      void setupFields();

      void updateDof(const char *name, const pylith::topology::FieldBase::DomainEnum domain, const int fiberDim);

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
      
      /** Add two fields, storing the result in one of the fields.
       *
       * @param field Field to add.
       */
      Field& operator+=(const Field& field);
      
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
       * Sieve section view.
       *
       * @param mesh Mesh associated with scatter.
       * @param context Label for context associated with vector.
       */
      template<typename scatter_mesh_type>
      void createScatter(const scatter_mesh_type& mesh,
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

      /// Scatter section information across processors to update the
      /// PETSc vector view of the field.
      void scatterSectionToVector(const char* context ="") const;

      /** Scatter section information across processors to update the
       * PETSc vector view of the field.
       *
       * @param vector PETSc vector to update.
       */
      void scatterSectionToVector(const PetscVec vector,
				  const char* context ="") const;

      /// Scatter PETSc vector information across processors to update the
      /// Sieve section view of the field.
      void scatterVectorToSection(const char* context ="") const;

      /** Scatter section information across processors to update the
       * PETSc vector view of the field.
       *
       * @param vector PETSc vector used in update.
       */
      void scatterVectorToSection(const PetscVec vector,
				  const char* context ="") const;

      /// Setup split field with all entries set to a default space of 0.
      void splitDefault(void);

    }; // Field

  } // topology
} // pylith


// End of file 
