// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/KinSrcLiuCos.i
 *
 * @brief Python interface to C++ KinSrcLiuCos object.
 */

namespace pylith {
    namespace faults {

	class KinSrcLiuCos : public pylith::faults::KinSrc {

	    // PUBLIC METHODS /////////////////////////////////////////////////
	public :

	    /// Default constructor.
	    KinSrcLiuCos(void);
      
	    /// Destructor.
	    ~KinSrcLiuCos(void);
      
	    // PROTECTED METHODS //////////////////////////////////////////////////
	protected:

	    /** Setup auxiliary subfields (discretization and query fns).
	     *
	     * @param[in] normalizer Normalizer for nondimensionalizing values.
	     * @param[in] cs Coordinate system for problem.
	     */
	    void _auxFieldSetup(const spatialdata::units::Nondimensional& normalizer,
				const spatialdata::geocoords::CoordSys* cs);
	    
	}; // class KinSrcLiuCos

    } // faults
} // pylith


// End of file 
