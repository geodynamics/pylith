// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file libsrc/meshio/OutputSubfield.hh
 *
 * @brief Manager for extracting subfield for output from a PyLith field.
 */

#if !defined(pylith_meshio_outputsubfield_hh)
#define pylith_meshio_outputsubfield_hh

#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/GenericComponent.hh" // ISA PyreComponent

#include "pylith/topology/FieldBase.hh" // HASA Description, Discretization

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/petscfwd.h" // HASA PetscVec

class pylith::meshio::OutputSubfield : public pylith::utils::GenericComponent {
    friend class TestOutputSubfield; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Create OutputSubfield from Field.
     *
     * @param[in] field Field with subfields.
     * @param[in] mesh Mesh for subfield.
     * @param[in] name Name of subfield that will be extracted.
     * @param[in] basisOrder Basis order for subfield.
     */
    static
    OutputSubfield* create(const pylith::topology::Field& field,
                           const pylith::topology::Mesh& mesh,
                           const char* name,
                           const int basisOrder);

    /** Create OutputSubfield from Field.
     *
     * @note Use this method in combination with extractSubfield() when there is no projection.
     *
     * @param[in] field Field with subfields.
     * @param[in] mesh Mesh for subfield.
     * @param[in] name Name of subfield that will be extracted.
     */
    static
    OutputSubfield* create(const pylith::topology::Field& field,
                           const pylith::topology::Mesh& mesh,
                           const char* name);

    /// Destructor
    ~OutputSubfield(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set label name and value.
     *
     * @param[in] name Name of PETSc label associated with subfield.
     * @param[in] value Value of PETSc label associated with subfield.
     */
    void setLabel(const char* name,
                  const int value);

    /** Get description of subfield.
     *
     * @returns Description of subfield.
     */
    const pylith::topology::FieldBase::Description& getDescription(void) const;

    /** Get basis order of subfield.
     *
     * @returns Basis order of subfield.
     */
    int getBasisOrder(void) const;

    /** Get PETSc global vector for projected subfield.
     *
     * @returns PETSc global vector.
     */
    PetscVec getVector(void) const;

    /** Get PETSc DM for projected subfield.
     *
     * @returns PETSc DM.
     */
    PetscDM getDM(void) const;

    /** Project PETSc vector to subfield.
     *
     * @param[in] fieldVector PETSc vector with subfields.
     */
    void project(const PetscVec& fieldVector);

    /** Project PETSc vector to subfield using label.
     *
     * @param[in] fieldVector PETSc vector with subfields.
     */
    void projectWithLabel(const PetscVec& fieldVector);

    /** Extract subfield from field.
     *
     * @pre DM must match for field and subfield.
     *
     * @param[in] field Field with all subfields.
     * @param[in] subfieldIndex Index of subfield to extract.
     */
    void extractSubfield(const pylith::topology::Field& field,
                         const PetscInt subfieldIndex);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    // Constructor.
    OutputSubfield(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::FieldBase::Description _description; ///< Description of subfield.
    pylith::topology::FieldBase::Discretization _discretization; ///< Discretization of subfield.
    PetscDM _dm; ///< PETSc DM for subfield.
    PetscVec _vector; ///< PETSc global vector for subfield.
    PetscPointFunc _fn; ///< PETSc point function for projection.
    PetscInt _subfieldIndex; ///< Index of subfield in fields.

    PetscDMLabel _label; ///< PETSc label associated with subfield.
    PetscInt _labelValue; ///< Value of PETSc label associated with subfield.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    OutputSubfield(const OutputSubfield&); ///< Not implemented.
    const OutputSubfield& operator=(const OutputSubfield&); ///< Not implemented

}; // OutputSubfield

#endif // pylith_meshio_outputsubfield_hh

// End of file
