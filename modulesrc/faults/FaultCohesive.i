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
        class FaultCohesive: public pylith::problems::Physics {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            FaultCohesive(void);

            /// Destructor.
            virtual ~FaultCohesive(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set name of label identifying cohesive cells.
             *
             * @param[in] value Name of label.
             */
            void setCohesiveLabelName(const char* value);

            /** Get name of label identifying cohesive cells.
             *
             * @returns Name of label.
             */
            const char* getCohesiveLabelName(void) const;

            /** Set identifier for fault cohesive cells.
             *
             * @param[in] value Fault identifier
             */
            void setCohesiveLabelValue(const int value);

            /** Get identifier for fault cohesive cells.
             *
             * @returns Fault identifier
             */
            int getCohesiveLabelValue(void) const;

            /** Set label marking surface of interface.
             *
             * @param[in] value Label of surface (from mesh generator).
             */
            void setSurfaceLabelName(const char* value);

            /** Get label marking surface of interface.
             *
             * @returns Label of surface (from mesh generator).
             */
            const char* getSurfaceLabelName(void) const;

            /** Set value of label marking surface of interface.
             *
             * @param[in] value Value of label of surface (from mesh generator).
             */
            void setSurfaceLabelValue(const int value);

            /** Get value of label marking surface of interface.
             *
             * @returns Value of label of surface (from mesh generator).
             */
            int getSurfaceLabelValue(void) const;

            /** Set label marking buried edges of interface surface.
             *
             * @param[in] value Label of buried surface edge (from mesh generator).
             */
            void setBuriedEdgesLabelName(const char* value);

            /** Get label marking buried edges of interface surface.
             *
             * @returns Label of buried surface edge (from mesh generator).
             */
            const char* getBuriedEdgesLabelName(void) const;

            /** Set value of label marking buried edges of interface surface.
             *
             * @param[in] value Value of label of buried surface edge (from mesh generator).
             */
            void setBuriedEdgesLabelValue(const int value);

            /** Get value of label marking buried edges of interface surface.
             *
             * @returns Value of label of buried surface edge (from mesh generator).
             */
            int getBuriedEdgesLabelValue(void) const;

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
