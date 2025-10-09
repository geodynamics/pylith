// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
/** @file modulesrc/faults/KinSrcRamp.i
 *
 * @brief Python interface to C++ KinSrcRamp object.
 */

namespace pylith {
    namespace faults {
        class KinSrcRamp: public pylith::faults::KinSrc {
            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Default constructor.
            KinSrcRamp(void);

            /// Destructor.
            ~KinSrcRamp(void);

            // PROTECTED METHODS //////////////////////////////////////////////////
protected:

            /** Setup auxiliary subfields (discretization and query fns).
             *
             * @param[in] scales Scales for nondimensionalizing values.
             * @param[in] cs Coordinate system for problem.
             */
            void _auxiliaryFieldSetup(const spatialdata::units::Scales& scales,
                                      const spatialdata::geocoords::CoordSys* cs);

        }; // class KinSrcRamp

    } // faults
} // pylith

// End of file
