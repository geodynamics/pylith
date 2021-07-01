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
 * @file modulesrc/topology/Mesh.hh
 *
 * @brief Python interface to C++ PyLith Mesh object.
 */

namespace pylith {
    namespace topology {
        // Mesh -------------------------------------------------------------
        class Mesh { // Mesh
                     // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Default constructor.
            Mesh(void);

            /** Constructor with dimension and communicator.
             *
             * @param dim Dimension associated with mesh cells.
             * @param comm MPI communicator for mesh.
             */
            Mesh(const int dim,
                 const MPI_Comm& comm=PETSC_COMM_WORLD);

            /// Default destructor
            ~Mesh(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Create clone.
             *
             * @returns Clone of mesh.
             */
            Mesh* clone(void) const;

            /** Get DMPlex mesh.
             *
             * @returns DMPlex mesh.
             */
            PetscDM getDM(void) const;

            /** Set DMPlex mesh.
             *
             * @param DMPlex mesh.
             * @param label Label for mesh.
             */
            void setDM(PetscDM dm,
                       const char* label="domain");

            /** Set coordinate system.
             *
             * @param cs Coordinate system.
             */
            void setCoordSys(const spatialdata::geocoords::CoordSys* cs);

            /** Get coordinate system.
             *
             * @returns Coordinate system.
             */
            const spatialdata::geocoords::CoordSys* getCoordSys(void) const;

            /** Get dimension of mesh.
             *
             * @returns Dimension of mesh.
             */
            int getDimension(void) const;

            /** Get MPI communicator associated with mesh.
             *
             * @returns MPI communicator.
             */
            const MPI_Comm getComm(void) const;

            /** Get MPI rank.
             *
             * @returns MPI rank.
             */
            int getCommRank(void) const;

            /** View mesh.
             *
             * @param viewOption PETSc DM view option.
             *
             * PETSc mesh view options include:
             *   short summary [empty]
             *   detail summary ::ascii_info_detail
             *   detail in a file :refined.mesh:ascii_info_detail
             *   latex in a file  :refined.tex:ascii_latex
             *   VTK vtk:refined.vtk:ascii_vtk
             */
            void view(const char* viewOption="") const;

            /** Get name of label for all mesh cells, including hybrid cells.
             *
             * @returns Name of label.
             */
            static
            const char* const getCellsLabelName(void);

        }; // Mesh

    } // topology
} // pylith

// End of file
