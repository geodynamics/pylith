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

/** @file modulesrc/faults/FaultCohesive.i
 *
 * @brief Python interface to C++ FaultCohesive object.
 */

namespace pylith {
    namespace faults {
        class FaultCohesive : public pylith::problems::Physics {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            FaultCohesive(void);

            /// Destructor.
            virtual ~FaultCohesive(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set identifier for fault cohesive cells.
             *
             * @param[in] value Fault identifier
             */
            void setInterfaceId(const int value);

            /** Get identifier for fault cohesive cells.
             *
             * @returns Fault identifier
             */
            int getInterfaceId(void) const;

            /** Set label marking surface of interface.
             *
             * @param[in] value Label of surface (from mesh generator).
             */
            void setSurfaceMarkerLabel(const char* value);

            /** Get label marking surface of interface.
             *
             * @returns Label of surface (from mesh generator).
             */
            const char* getSurfaceMarkerLabel(void) const;

            /** Set label marking buried edges of interface surface.
             *
             * @param[in] value Label of buried surface edge (from mesh generator).
             */
            void setBuriedEdgesMarkerLabel(const char* value);

            /** Get label marking buried edges of interface surface.
             *
             * @returns Label of buried surface edge (from mesh generator).
             */
            const char* getBuriedEdgesMarkerLabel(void) const;

            /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
             *
             * @param vec Reference direction unit vector.
             */
            void setRefDir1(const PylithReal vec[3]);

            /** Set second choice for reference direction to discriminate among tangential directions in 3-D.
             *
             * @param vec Reference direction unit vector.
             */
            void setRefDir2(const PylithReal vec[3]);

            /** Adjust mesh topology for fault implementation.
             *
             * @param mesh[in] PETSc mesh.
             */
            void adjustTopology(pylith::topology::Mesh* const mesh);

        }; // class FaultCohesive

    } // faults
} // pylith

// End of file
