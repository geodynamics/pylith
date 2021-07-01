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

/** @file modulesrc/bc/BoundaryCondition.i
 *
 * @brief Python interface to C++ BoundaryCondition object.
 */

namespace pylith {
    namespace bc {
        class BoundaryCondition : public pylith::problems::Physics {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            BoundaryCondition(void);

            /// Destructor.
            virtual ~BoundaryCondition(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set label marking boundary associated with boundary condition surface.
             *
             * @param[in] value Label of surface (from mesh generator).
             */
            void setMarkerLabel(const char* value);

            /** Get label marking boundary associated with boundary condition surface.
             *
             * @returns Label of surface (from mesh generator).
             */
            const char* getMarkerLabel(void) const;

            /** Set name of solution subfield associated with boundary condition.
             *
             * @param[in] value Name of solution subfield.
             */
            void setSubfieldName(const char* value);

            /** Get name of solution subfield associated with boundary condition.
             *
             * @preturn Name of solution subfield.
             */
            const char* getSubfieldName(void) const;

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

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void verifyConfiguration(const pylith::topology::Field& solution) const;

        };

        // class BoundaryCondition

    } // bc
} // pylith

// End of file
