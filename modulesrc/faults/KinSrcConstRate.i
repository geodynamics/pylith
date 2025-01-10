// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
/** @file modulesrc/faults/KinSrcConstRate.i
 *
 * @brief Python interface to C++ KinSrcConstRate object.
 */

namespace pylith {
    namespace faults {
        class KinSrcConstRate: public pylith::faults::KinSrc {
            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Default constructor.
            KinSrcConstRate(void);

            /// Destructor.
            ~KinSrcConstRate(void);

            // PROTECTED METHODS //////////////////////////////////////////////////
protected:

            /** Setup auxiliary subfields (discretization and query fns).
             *
             * @param[in] normalizer Normalizer for nondimensionalizing values.
             * @param[in] cs Coordinate system for problem.
             */
            void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                      const spatialdata::geocoords::CoordSys* cs);

        }; // class KinSrcConstRate

    } // faults
} // pylith

// End of file
