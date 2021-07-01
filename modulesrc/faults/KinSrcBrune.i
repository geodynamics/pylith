// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/KinSrcBrune.i
 *
 * @brief Python interface to C++ KinSrcBrune object.
 */

namespace pylith {
    namespace faults {

	class KinSrcBrune : public pylith::faults::KinSrc {

	    // PUBLIC METHODS /////////////////////////////////////////////////
	public :

	    /// Default constructor.
	    KinSrcBrune(void);
      
	    /// Destructor.
	    ~KinSrcBrune(void);
      
	    // PROTECTED METHODS //////////////////////////////////////////////////
	protected:

	    /** Setup auxiliary subfields (discretization and query fns).
	     *
	     * @param[in] normalizer Normalizer for nondimensionalizing values.
	     * @param[in] cs Coordinate system for problem.
	     */
	    void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
				const spatialdata::geocoords::CoordSys* cs);
	    
	}; // class KinSrcBrune

    } // faults
} // pylith


// End of file 
