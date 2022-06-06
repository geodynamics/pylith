// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

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
