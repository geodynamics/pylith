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
 * @file modulesrc/topology/Fields.i
 *
 * @brief Python interface to C++ Fields object.
 */

namespace pylith {
  namespace topology {

    template<typename field_type>
    class Fields
    { // Fields

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param mesh Finite-element mesh.
       */
      Fields(const typename field_type::Mesh& mesh);

      /// Destructor.
      virtual
      ~Fields(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);

      /** Check if fields contains a given field.
       *
       * @param name Name of field.
       * @return True if fields contains field, false otherwise.
       */
      bool hasField(const char* name) const;

      /** Add field.
       *
       * @param name Name of field.
       * @param label Label for field.
       */
      void add(const char* name,
	       const char* label);

      /** Add field.
       *
       * @param name Name of field.
       * @param label Label for field.
       * @param domain Type of points over which to define field.
       * @param fiberDim Fiber dimension for field.
       */
      void add(const char* name,
	       const char* label,
	       const pylith::topology::FieldBase::DomainEnum domain,
	       const int fiberDim);

      /** Delete field.
       *
       * @param name Name of field.
       */
      void delField(const char* name);

      /** Get field.
       *
       * @param name Name of field.
       */
      const field_type& get(const char* name) const;
	   
      /** Get field.
       *
       * @param name Name of field.
       */
      field_type& get(const char* name);
	   
      /** Copy layout to other fields.
       *
       * @param name Name of field to use as template for layout.
       */
      void copyLayout(const char* name);
      
      /** Get mesh associated with fields.
       *
       * @returns Finite-element mesh.
       */
      const typename field_type::Mesh& mesh(void) const;

      /** Return the names of all fields.
       *
       * @returns an array of all field names
       */
      void fieldNames(int *numNames, char ***outNames);

    }; // Fields

  } // topology
} // pylith


// End of file 
