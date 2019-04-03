// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
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

            /** Set label for marker associated with patch.
             *
             * @param[in] value Label for marker associated with patch.
             */
            void setMarkerLabel(const char* value);

            /** Get label for marker associated with patch.
             *
             * @returns Label for marker associated with patch.
             */
            const char* getMarkerLabel(void) const;

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
