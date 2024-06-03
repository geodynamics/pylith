// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/problems/Problem.hh
 *
 * @brief Python interface to C++ Problem.
 */

namespace pylith {
    namespace problems {
        class InitialConditionPatch: public pylith::problems::InitialCondition {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            InitialConditionPatch(void);

            /// Destructor
            virtual ~InitialConditionPatch(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set name of label marking material.
             *
             * @param[in] value Name of label for material (from mesh generator).
             */
            void setLabelName(const char* value);

            /** Get name of label marking material.
             *
             * @returns Name of label for material (from mesh generator).
             */
            const char* getLabelName(void) const;

            /** Set value of label marking material.
             *
             * @param[in] value Value of label for material (from mesh generator).
             */
            void setLabelValue(const int value);

            /** Get value of label marking material.
             *
             * @returns Value of label for material (from mesh generator).
             */
            int getLabelValue(void) const;

            /** Set spatial database holding initial conditions.
             *
             * @param[in] db Spatial database holding initial conditions.
             */
            void setDB(spatialdata::spatialdb::SpatialDB* db);

            /** Set solver type.
             *
             * @param[out] solution Solution field.
             * @param[in] normalizer Nondimensionalization.
             */
            void setValues(pylith::topology::Field* solution,
                           const spatialdata::units::Nondimensional& normalizer);

        }; // InitialConditionPatch

    } // problems
} // pylith

// End of file
