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

/** @file modulesrc/bc/TimeDependent.i
 *
 * @brief Python interface to C++ TimeDependent object.
 */

namespace pylith {
  namespace bc {

    class pylith::bc::TimeDependent
    { // class TimeDependent

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      TimeDependent(void);

      /// Destructor.
      ~TimeDependent(void);
      

      /** Set indices of degrees of freedom associated with BC.
       *
       * Note: Forces at all points are applied to the same degrees of freedom.
       *
       * Example: [0, 1] to apply forces to x and y degrees of freedom in
       * Cartesian system.
       *
       * @param flags Array of indices for degrees of freedom for forces.
       * @param size Size of array
       */
      %apply(int* INPLACE_ARRAY1, int DIM1) {
	(const int* flags, 
	 const int size)
	  };
      void bcDOF(const int* flags,
		 const int size);  
      %clear(const int* flags, const int size);
      
      /** Set database for initial values.
       *
       * @param db Spatial database
       */
      void dbInitial(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set database for rate of change of values.
       *
       * @param db Spatial database
       */
      void dbRate(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set database for change in values.
       *
       * @param db Spatial database
       */
      void dbChange(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set database for temporal evolution of change in value.
       *
       * @param db Time history database.
       */
      void dbTimeHistory(spatialdata::spatialdb::TimeHistory* const db);
      
      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      virtual
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;

      // PROTECTED METHODS //////////////////////////////////////////////////
    protected :
      
      /** Get label of boundary condition surface.
       *
       * @returns Label of surface (from mesh generator).
       */
      virtual
      const char* _getLabel(void) const = 0;
      
      /** Get manager of scales used to nondimensionalize problem.
       *
       * @returns Nondimensionalizer.
       */
      virtual
      const spatialdata::units::Nondimensional& _getNormalizer(void) const = 0;
      
    }; // class TimeDependent

  } // bc
} // pylith


// End of file 
