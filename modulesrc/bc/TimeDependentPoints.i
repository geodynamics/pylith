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

/** @file modulesrc/bc/TimeDependentPoints.i
 *
 * @brief Python interface to C++ TimeDependentPoints object.
 */

namespace pylith {
  namespace bc {

    class pylith::bc::TimeDependentPoints : public BoundaryConditionPoints, 
					    public TimeDependent
    { // class TimeDependentPoints

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :
      
      /// Default constructor.
      TimeDependentPoints(void);
      
      /// Destructor.
      ~TimeDependentPoints(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      // PROTECTED METHODS //////////////////////////////////////////////
    protected :
      
      /** Get label of boundary condition surface.
       *
       * @returns Label of surface (from mesh generator).
       */
      const char* _getLabel(void) const;

    }; // class TimeDependentPoints

  } // bc
} // pylith


// End of file 
