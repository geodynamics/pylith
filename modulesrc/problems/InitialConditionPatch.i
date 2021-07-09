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
// Copyright (c) 2010-2021 University of California, Davis
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
        class InitialConditionPatch : public pylith::problems::InitialCondition {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            InitialConditionPatch(void);

            /// Destructor
            virtual ~InitialConditionPatch(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set material id associated with patch.
             *
             * @param[in] value Material id associated with patch.
             */
            void setMaterialId(const int value);

            /** Get material id associated with patch.
             *
             * @returns Material id associated with patch.
             */
            int getMaterialId(void) const;

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
