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

    class Field : public FieldBase
    { // Field

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param mesh Finite-element mesh.
       */
      Field(const Mesh& mesh);

      /// Destructor.
      ~Field(void);
      
      /** Get mesh associated with field.
       *
       * @returns Finite-element mesh.
       */
      const Mesh& mesh(void) const;

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
      
    }; // Field
    
  } // topology
} // pylith


// End of file 
