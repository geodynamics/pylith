// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/materials/Material.i
 *
 * Python interface to C++ abstract base Material.
 */

namespace pylith {
  namespace materials {

    class Material
    { // class Material

    public:
      typedef pylith::topology::Field<pylith::topology::Mesh> field_type;
      // PUBLIC METHODS /////////////////////////////////////////////////
    public :
      
      /** Default constructor.
       *
       * @param dimension Spatial dimension associated with material.
       * @param tensorSize Array of names of database values for material.
       * @param metadata Metadata for physical properties and state variables.
       */
      Material(const int dimension,
	       const int tensorSize,
	       const pylith::materials::Metadata& metadata);

      /// Destructor.
      virtual
      ~Material(void);

      /** Get spatial dimension of material.
       *
       * @returns Spatial dimension.
       */
      int dimension(void) const;
      
      /** Set identifier of material.
       *
       * @param value Material identifier
       */
      void id(const int value);
      
      /** Get identifier of material.
       *
       * @returns Material identifier
       */
      int id(void) const;
      
      /** Set label of material.
       *
       * @param value Label of material
       */
      void label(const char* value);
      
      /** Get label of material.
       *
       * @returns Label of material
       */
      const char* label(void) const;
      
      /** Set current time step.
       *
       * @param dt Current time step.
       */
      virtual
      void timeStep(const double dt);
      
      /** Get current time step.
       *
       * @returns Current time step.
       */
      double timeStep(void) const;
      
      /** Set database for physical property parameters.
       *
       * @param value Pointer to database.
       */
      void dbProperties(spatialdata::spatialdb::SpatialDB* value);
      
      /** Set database for initial state variables.
       *
       * @param value Pointer to database.
       */
      void dbInitialState(spatialdata::spatialdb::SpatialDB* value);
      
      /** Set scales used to nondimensionalize physical properties.
       *
       * @param dim Nondimensionalizer
       */
      void normalizer(const spatialdata::units::Nondimensional& dim);
      
      /** Get size of stress/strain tensor associated with material.
       *
       * @returns Size of array holding stress/strain tensor.
       */
      int tensorSize(void) const;
      
      /** Get flag indicating whether Jacobian matrix must be reformed for
       * current state.
       *
       * @returns True if Jacobian matrix must be reformed, false otherwise.
       */
      bool needNewJacobian(void) const;
      
      /// Reset flag indicating whether Jacobian matrix must be reformed for
      /// current state.
      void resetNeedNewJacobian(void);
      
      /** Get physical property or state variable field. Data is returned
       * via the argument.
       *
       * @param field Field over material cells.
       * @param name Name of field to retrieve.
       */
      void getField(field_type *field, const char* name) const;

      /** Get the properties.
       *
       * @returns the properties
       */
      const field_type& getProperties() const;

      /** Get the state variables.
       *
       * @returns the state variables
       */
      const field_type& getStateVars() const;

      // PROTECTED METHODS //////////////////////////////////////////////
    protected :
      
      /** Compute properties from values in spatial database.
       *
       * @param propValues Array of property values.
       * @param dbValues Array of database values.
       */
      virtual
      void _dbToProperties(double* const propValues,
			   const double_array& dbValues) const = 0;
      
      /** Nondimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      virtual
      void _nondimProperties(double* const values,
			     const int nvalues) const = 0;
      
      /** Dimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      virtual
      void _dimProperties(double* const values,
			  const int nvalues) const = 0;

    }; // class Material

  } // materials
} // pylith


// End of file 
