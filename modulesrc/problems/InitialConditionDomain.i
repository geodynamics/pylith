// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


namespace pylith {
    namespace problems {
        class InitialConditionDomain : public pylith::problems::InitialCondition {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
          public:
            /// Constructor
            InitialConditionDomain(void);

            /// Destructor
            virtual ~InitialConditionDomain(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set spatial database holding initial conditions.
             *
             * @param[in] db Spatial database holding initial conditions.
             */
            void setDB(spatialdata::spatialdb::SpatialDB* db);

            /** Set solver type.
             *
             * @param[out] solution Solution field.
             * @param[in] scales Nondimensionalization.
             */
            void setValues(pylith::topology::Field* solution,
                           const pylith::scales::Scales& scales);

        }; // InitialConditionDomain

    } // namespace problems
} // namespace pylith


// End of file
