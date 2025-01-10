// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/meshio/MeshIO.hh" // ISA MeshIO

#include <string> // HASA std::string

class pylith::meshio::MeshIOCubit : public MeshIO {
    friend class TestMeshIOCubit; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIOCubit(void);

    /// Destructor
    ~MeshIOCubit(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set filename for Cubit file.
     *
     * @param filename Name of file
     */
    void setFilename(const char* name);

    /** Get filename of Cubit file.
     *
     * @returns Name of file
     */
    const char* getFilename(void) const;

    /** Set flag on whether to use nodeset ids or names.
     *
     * @param flag True to use node set names.
     */
    void setUseNodesetNames(const bool flag);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    void _write(void) const;

    /// Read mesh
    void _read(void);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Read mesh vertices.
     *
     * @param ncfile Cubit Exodus file.
     * @param coordinates Pointer to array of vertex coordinates.
     * @param numVertices Pointer to number of vertices.
     * @param spaceDim Pointer to dimension of coordinates vector space.
     */
    void _readVertices(ExodusII& filein,
                       scalar_array* coordinates,
                       int* numVertices,
                       int* spaceDim) const;

    /** Read mesh cells.
     *
     * @param ncfile Cubit Exodus file.
     * @param pCells Pointer to array of indices of cell vertices
     * @param pMaterialIds Pointer to array of material identifiers
     * @param pNumCells Pointer to number of cells
     * @param pNumCorners Pointer to number of corners
     */
    void _readCells(ExodusII& filein,
                    int_array* pCells,
                    int_array* pMaterialIds,
                    int* numCells,
                    int* numCorners) const;

    /** Read point groups.
     *
     * @param ncfile Cubit Exodus file.
     */
    void _readGroups(ExodusII& filein);

    /** Write mesh dimensions.
     *
     * @param ncfile Cubit Exodus file.
     */
    void _writeDimensions(ExodusII& ncfile) const;

    /** Write mesh variables.
     *
     * @param ncfile Cubit Exodus file.
     */
    void _writeVariables(ExodusII& ncfile) const;

    /** Write mesh attributes.
     *
     * @param ncfile Cubit Exodus file.
     */
    void _writeAttributes(ExodusII& ncfile) const;

    /** Reorder vertices in cells to match PyLith conventions.
     *
     * @param cells Array of vertex indices for each cell [numCells*numCorners].
     * @param numCells Number of cells.
     * @param numCorners Number of vertices per cell.
     * @param meshDim Spatial dimension of mesh.
     */
    static
    void _orientCells(int_array* const cells,
                      const int numCells,
                      const int numCorners,
                      const int meshDim);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::string _filename; ///< Name of file
    bool _useNodesetNames; ///< True to use node set names instead of ids.

}; // MeshIOCubit

#include "MeshIOCubit.icc" // inline methods

// End of file
