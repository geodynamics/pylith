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
 * @file modulesrc/topology/Field.hh
 *
 * @brief Python interface to C++ Field object.
 */

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
      void scale(const double value);

      /** Get scale for dimensionalizing field.
       *
       * @returns Scale associated with field.
       */
      double scale(void) const;
      
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
      
      /// Create sieve section.
      void newSection(void);

      /** Create sieve section and set chart and fiber dimesion.
       *
       * @param domain Type of points over which to define section.
       * @param dim Fiber dimension for section.
       */
      void newSection(const pylith::topology::FieldBase::DomainEnum domain,
		      const int fiberDim);

      /** Create section with same layout (fiber dimension and
       * constraints) as another section. This allows the layout data
       * structures to be reused across multiple fields, reducing memory
       * usage.
       *
       * @param sec Section defining layout.
       */
      void newSection(const Field& src);

      /// Clear variables associated with section.
      void clear(void);

      /// Allocate field.
      void allocate(void);
      
      /// Zero section values.
      void zero(void);
      
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
      void operator+=(const Field& field);
      
      /** Dimensionalize field. Throws runtime_error if field is not
       * allowed to be dimensionalized.
       */
      void dimensionalize(void);
      
      /** Print field to standard out.
       *
       * @param label Label for output.
       */
      void view(const char* label);

      /// Create PETSc vector for field.
      void createVector(void);
      
      /** Get PETSc vector associated with field.
       *
       * @returns PETSc vector.
       */
      PetscVec vector(void);
      
      /** Get PETSc vector associated with field.
       *
       * @returns PETSc vector.
       */
      const PetscVec vector(void) const;
      
      /// Create PETSc vector scatter for field. This is used to transfer
      /// information from the "global" PETSc vector view to the "local"
      /// Sieve section view.
      void createScatter(void);
      
      /// Scatter section information across processors to update the
      /// PETSc vector view of the field.
      void scatterSectionToVector(void) const;
      
      /// Scatter PETSc vector information across processors to update the
      /// Sieve section view of the field.
      void scatterVectorToSection(void) const;

    }; // Field

  } // topology
} // pylith


// End of file 
