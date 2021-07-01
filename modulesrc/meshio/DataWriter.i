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
 * @file modulesrc/meshio/DataWriter.i
 *
 * @brief Python interface to C++ DataWriter object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::DataWriter {
            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Constructor
            DataWriter(void);

            /// Destructor
            virtual ~DataWriter(void);

            /** Make copy of this object.
             *
             * @returns Copy of this.
             */
            virtual
            DataWriter* clone(void) const = 0;

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set time scale for simulation time.
             *
             * @param value Time scale
             */
            void setTimeScale(const PylithScalar value);

            /** Is data writer open, i.e., ready for openTimeStep()/closeTimeStep()?
             *
             * @returns True if data writer is open, false otherwise.
             */
            bool isOpen(void) const;

            /** Prepare for writing files.
             *
             * @param mesh Finite-element mesh.
             * @param isInfo True if only writing info values.
             */
            virtual
            void open(const pylith::topology::Mesh& mesh,
                      const bool isInfo);

            /// Close output files.
            virtual
            void close(void);

            /** Prepare file for data at a new time step.
             *
             * @param t Time stamp for new data
             * @param mesh PETSc mesh object
             */
            virtual
            void openTimeStep(const PylithScalar t,
                              const pylith::topology::Mesh& mesh);

            /// Cleanup after writing data for a time step.
            virtual
            void closeTimeStep(void);

            /** Write field over vertices to file.
             *
             * @param[in] t Time associated with field.
             * @param[in] subfield Subfield with basis order 1.
             */
            virtual
            void writeVertexField(const PylithScalar t,
                                  const pylith::meshio::OutputSubfield& field) = 0;

            /** Write field over cells to file.
             *
             * @param[in] t Time associated with field.
             * @param[in] subfield Subfield with basis order 0.
             */
            virtual
            void writeCellField(const PylithScalar t,
                                const pylith::meshio::OutputSubfield& subfield) = 0;

            /** Write dataset with names of points to file.
             *
             * @param names Array with name for each point, e.g., station name.
             * @param mesh Finite-element mesh.
             *
             * Primarily used with OutputSolnPoints.
             */
            virtual
            void writePointNames(const pylith::string_vector& names,
                                 const pylith::topology::Mesh& mesh);

        }; // DataWriter

    } // meshio
} // pylith

// End of file
