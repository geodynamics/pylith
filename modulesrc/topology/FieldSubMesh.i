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
 * @file modulesrc/topology/FieldSubMesh.hh
 *
 * @brief Python interface to C++ FieldSubMesh object.
 */

namespace pylith {
  namespace topology {

    class FieldSubMesh : public FieldBase
    { // FieldSubMesh

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param mesh Lower dimension finite-element mesh.
       */
      FieldSubMesh(const SubMesh& mesh);

      /// Destructor.
      ~FieldSubMesh(void);

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
      void newSection(const DomainEnum domain,
		      const int fiberDim);

      /** Create section with same layout (fiber dimension and
       * constraints) as another section. This allows the layout data
       * structures to be reused across multiple fields, reducing memory
       * usage.
       *
       * @param sec Section defining layout.
       */
      void newSection(const FieldSubMesh& src);
      
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
       * @param field FieldSubMesh to copy.
       */
      void copy(const FieldSubMesh& field);
      
      /** Add two fields, storing the result in one of the fields.
       *
       * @param field FieldSubMesh to add.
       */
      void operator+=(const FieldSubMesh& field);
      
      /** Dimensionalize field. Throws runtime_error if field is not
       * allowed to be dimensionalized.
       */
      void dimensionalize(void);
      
      /** Print field to standard out.
       *
       * @param label Label for output.
       */
      void view(const char* label);
      
    }; // FieldSubMesh

  } // topology
} // pylith


// End of file 
