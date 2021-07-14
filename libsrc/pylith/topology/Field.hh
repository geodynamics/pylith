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
 * @file libsrc/topology/Field.hh
 *
 * @brief Vector field over the vertices or cells of a finite-element
 * mesh.
 */

#if !defined(pylith_topology_field_hh)
#define pylith_topology_field_hh

#include "FieldBase.hh" // ISA FieldBase
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/utils/petscfwd.h" // HASA PetscVec
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys

#include <map> // USES std::map
#include <string> // USES std::string

/** @brief Vector field over the vertices or cells of a finite-element
 * mesh.
 *
 * Extends PETSc section and vector by adding metadata.
 *
 * General steps for setting up a field:
 *
 * 1. Add subfields (description and discretization).
 * 2. Set the chart.
 * 3. Set any constraints (known DOF).
 * 4. Allocate.
 *
 * For local fields associated with parameters, we do not have step 3.
 * Step 3 is handled automatically by the PETSc DS if
 * DMAddBoundary is called.
 */
class pylith::topology::Field : public pylith::topology::FieldBase, public pylith::utils::GenericComponent {
    friend class FieldQuery; // Fill field using data.

    friend class TestFieldMesh; // unit testing
    friend class TestFieldSubmesh; // unit testing

    // PUBLIC ENUMS ///////////////////////////////////////////////////////
public:

    enum ViewOptions {
        VIEW_METADATA=0, ///< View metadata only.
        VIEW_LAYOUT=1, ///< View metadata and section.
        VIEW_VALUES=2, ///< View metadata and vector.
        VIEW_ALL=3, ///< View metadata, section, and vector.
    };

    // PUBLIC STRUCTS /////////////////////////////////////////////////////////////////////////////
public:

    /// Subfield auxiliary information.
    struct SubfieldInfo {
        Description description; ///< Description for subfield.
        Discretization fe; ///< Discretization information for subfield.
        int index; ///< Index of subfield in field.
    }; // SubfieldInfo

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param mesh Finite-element mesh for field (will be cloned).
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
     * @returns Finite-element mesh associated with field.
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

    /** Get spatial dimension of coordinate system for field.
     *
     * @returns Spatial dimension.
     */
    size_t getSpaceDim(void) const;

    /** Get the number of points in the chart.
     *
     * @returns the chart size.
     */
    PylithInt getChartSize(void) const;

    /** Get the number of degrees of freedom.
     *
     * @returns the number of degrees of freedom.
     */
    PylithInt getStorageSize(void) const;

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

    /** Add subfield to current field.
     *
     * Should be followed by calls to subfieldsSetup() and allocate().
     *
     * @param[in] description Description for subfield.
     * @param[in] discretization Deiscretization information for subfield.
     */
    void subfieldAdd(const Description& description,
                     const Discretization& discretization);

    /** Setup sections for subfields.
     *
     * Should be preceded by calls to subfieldAdd() and followed by calls to allocate().
     *
     * Additional setup of the field, such as customization of the discretization, can be done after calling this
     * routine.
     */
    void subfieldsSetup(void);

    /** Does field have given subfield?
     *
     * @param name Name of subfield.
     * @returns True if field has given subfield.
     */
    bool hasSubfield(const char* name) const;

    /** Get names of subfields.
     *
     * @returns Array of names of subfields.
     */
    pylith::string_vector getSubfieldNames(void) const;

    /** Get auxiliary information for subfield.
     *
     * @param name Name of field.
     * @returns Auxiliary information (including metadata) for subfield.
     */
    const SubfieldInfo& getSubfieldInfo(const char* name) const;

    /** Display field.
     *
     * @param[in] label Label for output.
     * @param[in] options Viewing options.
     */
    void view(const char* label,
              const ViewOptions options=VIEW_ALL) const;

    /** Scatter section information across processors to update the
     * global view of the field.
     *
     * @param[out] vector PETSc vector to update.
     * @param[in] mode Mode for scatter (INSERT_VALUES, ADD_VALUES).
     */
    void scatterLocalToVector(const PetscVec vector,
                              InsertMode mode=INSERT_VALUES) const;

    /** Scatter global information across processors to update the local
     * view of the field.
     *
     * @param[in] vector PETSc vector used in update.
     * @param[in] mode Mode for scatter (INSERT_VALUES, ADD_VALUES).
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

    // PRIVATE TYPEDEFS ///////////////////////////////////////////////////////////////////////////
private:

    typedef std::map<std::string, SubfieldInfo> subfields_type;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    subfields_type _subfields; ///< Map of subfields in field.
    std::string _label; ///< Label for field.

    pylith::topology::Mesh* _mesh; ///< Mesh associated with field.
    PetscVec _localVec; ///< Local PETSc vector.
    PetscVec _globalVec; ///< Global PETSc vector.
    PetscVec _outputVec; ///< Global PETSc vector without constrained DOF for output.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    const Field& operator=(const Field&); ///< Not implemented

}; // Field

#endif // pylith_topology_field_hh

// End of file
