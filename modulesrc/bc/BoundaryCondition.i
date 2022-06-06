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
// Copyright (c) 2010-2022 University of California, Davis
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
        class BoundaryCondition: public pylith::problems::Physics {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            BoundaryCondition(void);

            /// Destructor.
            virtual ~BoundaryCondition(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

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

            /** Set name of label marking boundary associated with boundary condition surface.
             *
             * @param[in] value Name of label for surface (from mesh generator).
             */
            void setLabelName(const char* value);

            /** Get name of label marking boundary associated with boundary condition surface.
             *
             * @returns Name of label for surface (from mesh generator).
             */
            const char* getLabelName(void) const;

            /** Set value of label marking boundary associated with boundary condition surface.
             *
             * @param[in] value Value of label for surface (from mesh generator).
             */
            void setLabelValue(const int value);

            /** Get value of label marking boundary associated with boundary condition surface.
             *
             * @returns Value of label for surface (from mesh generator).
             */
            int getLabelValue(void) const;

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
