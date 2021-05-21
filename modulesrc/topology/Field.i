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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
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
            const pylith::topology::Mesh& mesh(void) const;

            /** Set label for field.
             *
             * @param value Label for field.
             */
            void setLabel(const char* value);

            /** Get label for field.
             *
             * @returns Label for field.
             */
            const char* getLabel(void) const;

            /** Get spatial dimension of domain.
             *
             * @returns Spatial dimension of domain.
             */
            int getSpaceDim(void) const;

            /** Get the number of points in the chart.
             *
             * @returns the chart size.
             */
            PetscInt chartSize(void) const;

            /** Get the number of degrees of freedom.
             *
             * @returns the number of degrees of freedom.
             */
            PetscInt getStorageSize(void) const;

            /** Add subfield to current field.
             *
             * Should be followed by calls to subfieldsSetup() and allocate().
             *
             * @param[in] name Programatic name for subfield.
             * @param[in] alias User-specified name for subfield.
             * @param[in] fieldType Type of vector field.
             * @param[in] components Names of components in subfield.
             * @param[in] numComponents Number of components in subfield.
             * @param[in] basisOrder Polynomial order for basis.
             * @param[in] quadOrder Order of quadrature rule.
             * @param[in] dimension Dimension of points for discretization.
             * @param[in] cellBasis Type of basis functions to use (e.g., simplex, tensor, or default).
             * @param[in] isBasisContinuous True if basis is continuous.
             * @param[in] feSpace Finite-element space (polynomial or point).
             * @param[in] scale Scale for dimensionalizing field.
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
                             const CellBasis cellBasis,
                             const bool isBasisContinuous,
                             const SpaceEnum feSpace);

            %clear(const char* components[], const int numComponents);

            /** Setup sections for subfields.
             *
             * Should be preceded by calls to subfieldAdd() and followed by calls to subfieldSetDof().
             */
            void subfieldsSetup(void);

            /// Allocate field.
            void allocate(void);

            /// Zero section values (including constrained DOF).
            void zeroLocal(void);

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
