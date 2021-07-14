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
 * @file modulesrc/topology/Field.hh
 *
 * @brief Python interface to C++ Field object.
 */

namespace pylith {
    namespace topology {
        class Field : public FieldBase {
            // PUBLIC ENUMS ///////////////////////////////////////////////////////
public:

            enum ViewOptions {
                VIEW_METADATA=0, ///< View metadata only.
                VIEW_LAYOUT=1, ///< View metadata and section.
                VIEW_VALUES=2, ///< View metadata and vector.
                VIEW_ALL=3, ///< View metadata, section, and vector.
            };

            // PUBLIC MEMBERS /////////////////////////////////////////////////
public:

            /** Default constructor.
             *
             * @param mesh Finite-element mesh.
             */
            Field(const pylith::topology::Mesh& mesh);

            /** Constructor with field to use for layout.
             *
             * @param[in] src Field to copy layout from.
             *
             * @note Don't forget to call setLabel(), especially if reusing a field.
             */
            Field(const Field& src);

            /// Destructor.
            ~Field(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get mesh associated with field.
             *
             * @returns Finite-element mesh.
             */
            const pylith::topology::Mesh& getMesh(void) const;

            /** Get PETSc DM associated with field.
             *
             * @returns PETSc DM
             */
            PetscDM getDM(void) const;

            /** Get label for field.
             *
             * @returns Label for field.
             */
            const char* getLabel(void) const;

            /** Set label for field.
             *
             * @param value Label for field.
             */
            void setLabel(const char* value);

            /** Get local PetscSection.
             *
             * @returns PETSc section.
             */
            PetscSection getLocalSection(void) const;

            /** Get global PetscSection.
             *
             * @returns PETSc section.
             */
            PetscSection getGlobalSection(void) const;

            /** Get the local PETSc Vec.
             *
             * @returns PETSc Vec object.
             */
            PetscVec getLocalVector(void) const;

            /** Get the global PETSc Vec.
             *
             * @returns PETSc Vec object.
             */
            PetscVec getGlobalVector(void) const;

            /** Get the global PETSc Vec without constrained degrees of freedom for output.
             *
             * @returns PETSc Vec object.
             */
            PetscVec getOutputVector(void) const;

            /** Get spatial dimension of domain.
             *
             * @returns Spatial dimension of domain.
             */
            int getSpaceDim(void) const;

            /** Get the number of points in the chart.
             *
             * @returns the chart size.
             */
            PetscInt getChartSize(void) const;

            /** Get the number of degrees of freedom.
             *
             * @returns the number of degrees of freedom.
             */
            PetscInt getStorageSize(void) const;

            /** Create discretization for field.
             *
             * @important Should be called for all fields after
             * Field::subfieldsSetup() and before PetscDSAddBoundary() and
             * Field::allocate().
             */
            void createDiscretization(void);

            /// Allocate field and zero the local vector.
            void allocate(void);

            /// Zero local values (including constrained values).
            void zeroLocal(void);

            /** Add subfield to current field (for use from SWIG).
             *
             * Should be followed by calls to subfieldsSetup() and allocate().
             *
             * @param[in] name Programatic name of subfield.
             * @param[in] alias User-specified alias for subfield.
             * @param[in] fieldType Type of vector field.
             * @param[in] components Array of names of field components.
             * @param[in] numComponents Size of array.
             * @param[in] scale Dimensional scale associated with field.
             * @param[in] basisOrder Order of basis functions for discretization.
             * @param[in] quadOrder Order of numerical quadrature for discretization.
             * @param[in] dimension Dimension of points for discretization.
             * @param[in] isFaultOnly True if subfield is limited to fault degrees of freedom.
             * @param[in] cellBasis Type of basis functions to use (e.g., simplex, tensor, or default).
             * @param[in] feSpace Finite-element space (POLYNOMIAL_SPACE or POINT_SPACE).
             * @param[in] isBasisContinuous True if basis is continuous.
             */
            %apply(const char* const* string_list, const int list_len) {
                (const char* components[], const int numComponents)
            };
            void subfieldAdd(const char *name,
                             const char* alias,
                             const VectorFieldEnum fieldType,
                             const char* components[],
                             const int numComponents,
                             const double scale,
                             const int basisOrder,
                             const int quadOrder,
                             const int dimension,
                             const bool isFaultOnly,
                             const CellBasis cellBasis,
                             const SpaceEnum feSpace,
                             const bool isBasisContinuous);

            %clear(const char* components[], const int numComponents);

            /** Setup sections for subfields.
             *
             * Should be preceded by calls to subfieldAdd() and followed by calls to subfieldSetDof().
             */
            void subfieldsSetup(void);

            /** Does field have given subfield?
             *
             * @param name Name of subfield.
             * @returns True if field has given subfield.
             */
            bool hasSubfield(const char* name) const;

            /** Print field to standard out.
             *
             * @param label Label for output.
             */
            void view(const char* label,
                      const ViewOptions options=VIEW_ALL);

            /** Scatter section information across processors to update the
             * global view of the field.
             *
             * @param vector PETSc vector to update.
             * @param context Label for context associated with vector.
             */
            void scatterLocalToVector(const PetscVec vector,
                                      InsertMode mode=INSERT_VALUES) const;

            /** Scatter global information across processors to update the local
             * view of the field.
             *
             * @param vector PETSc vector used in update.
             * @param context Label for context associated with vector.
             */
            void scatterVectorToLocal(const PetscVec vector,
                                      InsertMode mode=INSERT_VALUES) const;

            /** Scatter section information across processors to update the
             * output view of the field.
             *
             * @param[in] mode Mode for scatter (INSERT_VALUES, ADD_VALUES).
             */
            void scatterLocalToOutput(InsertMode mode=INSERT_VALUES) const;

            /// Create global vector.
            void createGlobalVector(void);

            /// Create global vector with no constrained degrees of freedom for output.
            void createOutputVector(void);

        }; // Field

    } // topology
} // pylith

// End of file
