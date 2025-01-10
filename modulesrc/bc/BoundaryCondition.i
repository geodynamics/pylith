// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
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

            /** Create diagnostic field.
             *
             * @param[in] solution Solution field.
             * @param[in] physicsMesh Finite-element mesh associated with physics.
             *
             * @returns Diagnostic field if applicable, otherwise NULL.
             */
            virtual
            pylith::topology::Field* createDiagnosticField(const pylith::topology::Field& solution,
                                                           const pylith::topology::Mesh& physicsMesh);

        };

        // class BoundaryCondition

    } // bc
} // pylith

// End of file
