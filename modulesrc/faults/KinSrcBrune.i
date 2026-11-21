// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


namespace pylith {
    namespace faults {
        class KinSrcBrune : public pylith::faults::KinSrc {
            // PUBLIC METHODS /////////////////////////////////////////////////
          public:
            /// Default constructor.
            KinSrcBrune(void);

            /// Destructor.
            ~KinSrcBrune(void);

            // PROTECTED METHODS //////////////////////////////////////////////////
          protected:
            /** Setup auxiliary subfields (discretization and query fns).
             *
             * @param[in] scales Scales for nondimensionalizing values.
             * @param[in] cs Coordinate system for problem.
             */
            void _auxiliaryFieldSetup(const pylith::scales::Scales& scales,
                                      const spatialdata::geocoords::CoordSys* cs);

        }; // class KinSrcBrune

    } // namespace faults
} // namespace pylith


// End of file
