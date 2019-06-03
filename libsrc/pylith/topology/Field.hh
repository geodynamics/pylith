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
 * @file libsrc/topology/Field.hh
 *
 * @brief Vector field over the vertices or cells of a finite-element
 * mesh.
 */

#if !defined(pylith_topology_field_hh)
#define pylith_topology_field_hh

// Include directives ---------------------------------------------------
#include "FieldBase.hh" // ISA FieldBase
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/utils/arrayfwd.hh" // HASA int_array
#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include <map> // USES std::map
#include <string> // USES std::string

// Field ----------------------------------------------------------------
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

    // PUBLIC STRUCTS ///////////////////////////////////////////////////////
public:

    /// Subfield auxiliary information.
    struct SubfieldInfo {
        Description description; ///< Description for subfield.
        Discretization fe; ///< Discretization information for subfield.
        int index; ///< Index of subfield in field.
        PetscDM dm; ///< PETSc DM associated with subfield.
    }; // SubfieldInfo

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param mesh Finite-element mesh.
     */
    Field(const Mesh& mesh);

    /** Constructor with mesh, PETSc DM, local data, and metadata.
     *
     * @param mesh Finite-element mesh.
     * @param dm PETSc DM for field.
     * @param localVec PETSc Vec with local data for field.
     * @param description Field description.
     */
    Field(const Mesh& mesh,
          PetscDM dm,
          PetscVec localVec,
          const Description& description);

    /// Destructor.
    ~Field(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get mesh associated with field.
     *
     * @returns Finite-element mesh.
     */
    const Mesh& mesh(void) const;

    /** Get PETSc DM associated with field.
     *
     * @returns PETSc DM
     */
    PetscDM dmMesh(void) const;

    /** Set label for field.
     *
     * @param value Label for field.
     */
    void label(const char* value);

    /** Get label for field.
     *
     * @returns Label for field.
     */
    const char* label(void) const;

    /** Set flag indicating whether it is okay to dimensionalize field.
     *
     * @param value True if it is okay to dimensionalize field.
     */
    void dimensionalizeOkay(const bool value);

    /** Set flag indicating whether it is okay to dimensionalize field.
     *
     * @param value True if it is okay to dimensionalize field.
     */
    bool dimensionalizeOkay(void) const;

    /** Get spatial dimension of domain.
     *
     * @returns Spatial dimension of domain.
     */
    int spaceDim(void) const;

    /** Get vector field type.
     *
     * If field constains only one subfield, return the vector
     * field type for that subfield, otherwise return OTHER.
     */
    VectorFieldEnum vectorFieldType(void) const;

    /** Get the number of points in the chart.
     *
     * @returns the chart size.
     */
    PylithInt chartSize(void) const;

    /** Get the number of degrees of freedom.
     *
     * @returns the number of degrees of freedom.
     */
    PylithInt sectionSize(void) const;

    /** Get local PetscSection.
     *
     * @returns PETSc section.
     */
    PetscSection localSection(void) const;

    /** Get global PetscSection.
     *
     * @returns PETSc section.
     */
    PetscSection globalSection(void) const;

    /** Get the local PETSc Vec.
     *
     * @returns PETSc Vec object.
     */
    PetscVec localVector(void) const;

    /** Create section with same layout (fiber dimension and
     * constraints) as another section. This allows the layout data
     * structures to be reused across multiple fields, reducing memory
     * usage.
     *
     * @param src Field defining layout.
     *
     * @note Don't forget to call label(), especially if reusing a field.
     */
    void cloneSection(const Field& src);

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
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space (POLYNOMIAL_SPACE or POINT_SPACE).
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
                     const bool isBasisContinuous,
                     const SpaceEnum feSpace);

    /** Add subfield to current field.
     *
     * Should be followed by calls to subfieldsSetup() and allocate().
     *
     * @param[in] description Description for subfield.
     * @param[in] discretization Deiscretization information for subfield.
     */
    void subfieldAdd(const Description& description,
                     const Discretization& discretization);

    /** Update subfield description and discretization.
     *
     * @param[in] subfieldName Name of subfield to update.
     * @param[in] subfieldNameNew New name for subfield.
     * @param[in] description Description for new subfield.
     * @param[in] discretization Deiscretization information for new subfield.
     */
    void subfieldUpdate(const char* subfieldName,
                        const char* subfieldNameNew,
                        const Description& description,
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
    pylith::string_vector subfieldNames(void) const;

    /** Get auxiliary information for subfield.
     *
     * @param name Name of field.
     * @returns Auxiliary information (including metadata) for subfield.
     */
    const SubfieldInfo& subfieldInfo(const char* name) const;

    /// Clear variables associated with section.
    void clear(void);

    /// Allocate field and zero local vector..
    void allocate(void);

    /// Zero local values (including constrained values).
    void zeroLocal(void);

    /** Copy field values and metadata.
     *
     * @param field Field to copy.
     */
    void copy(const Field& field);

    /** Copy subfield values and its metadata to field;
     *
     * @param field Field to copy from.
     * @param name Name of subfield to copy.
     */
    void copySubfield(const Field& field,
                      const char* name);

    /** Dimensionalize field. Throws runtime_error if field is not
     * allowed to be dimensionalized.
     */
    void dimensionalize(void) const;

    /** Print field to standard out.
     *
     * @param[in] label Label for output.
     * @param[in] options Viewing options.
     */
    void view(const char* label,
              const ViewOptions options=VIEW_ALL) const;

    /** Create PETSc vector scatter for field. This is used to transfer
     * information from the "global" PETSc vector view to the "local"
     * PETSc section view. The PETSc vector does not contain constrained
     * DOF. Use createScatterWithBC() to include the constrained DOF in
     * the PETSc vector.
     *
     * @param mesh Mesh associated with scatter.
     * @param context Label for context associated with vector.
     */
    void createScatter(const Mesh& mesh,
                       const char* context);

    /** Create PETSc vector scatter for field. This is used to transfer
     * information from the "global" PETSc vector view to the "local"
     * PETSc section view. The PETSc vector includes constrained
     * DOF. Use createScatter() if constrained DOF should be omitted
     * from the PETSc vector.
     *
     * @param context Label for context associated with vector.
     */
    void createScatterWithBC(const char* context);

    /** Create PETSc vector scatter for field. This is used to transfer
     * information from the "global" PETSc vector view to the "local"
     * PETSc section view. The PETSc vector includes constrained
     * DOF. Use createScatter() if constrained DOF should be omitted
     * from the PETSc vector.
     *
     * @param mesh Mesh associated with scatter.
     * @param context Label for context associated with vector.
     */
    void createScatterWithBC(const Mesh& mesh,
                             const char* context);

    /** Get PETSc vector associated with scatter for field.
     *
     * @param context Label for context associated with vector.
     * @returns PETSc vector.
     */
    PetscVec scatterVector(const char* context);

    /** Get PETSc vector associated with scatter for field.
     *
     * @param context Label for context associated with vector.
     * @returns PETSc vector.
     */
    const PetscVec scatterVector(const char* context) const;

    /** Scatter section information across processors to update the
     * global view of the field.
     *
     * @param[in] context Label for context associated with vector.
     * @param[in] mode Mode for scatter (INSERT_VALUES, ADD_VALUES).
     */
    void scatterLocalToContext(const char* context,
                               InsertMode mode=INSERT_VALUES) const;

    /** Scatter section information across processors to update the
     * global view of the field.
     *
     * @param[out] vector PETSc vector to update.
     * @param[in] context Label for context associated with vector.
     * @param[in] mode Mode for scatter (INSERT_VALUES, ADD_VALUES).
     */
    void scatterLocalToVector(const PetscVec vector,
                              InsertMode mode=INSERT_VALUES) const;

    /** Scatter global information across processors to update the local
     * view of the field.
     *
     * @param[in] context Label for context associated with vector.
     * @param[in] mode Mode for scatter (INSERT_VALUES, ADD_VALUES).
     */
    void scatterContextToLocal(const char* context,
                               InsertMode mode=INSERT_VALUES) const;

    /** Scatter global information across processors to update the local
     * view of the field.
     *
     * @param[in] vector PETSc vector used in update.
     * @param[in] context Label for context associated with vector.
     * @param[in] mode Mode for scatter (INSERT_VALUES, ADD_VALUES).
     */
    void scatterVectorToLocal(const PetscVec vector,
                              InsertMode mode=INSERT_VALUES) const;

    // PRIVATE STRUCTS //////////////////////////////////////////////////////
private:

    /// Data structures used in scattering to/from PETSc Vecs.
    struct ScatterInfo {
        PetscDM dm; ///< PETSc DM defining the communication pattern
        PetscVec vector; ///< PETSc vector associated with field.
    }; // ScatterInfo

    // PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private:

    typedef std::map<std::string, ScatterInfo> scatter_map_type;
    typedef std::map<std::string, SubfieldInfo> subfields_type;

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    /** Setup field to hold values extract from subfield of another field.
     *
     * @param field Field to copy from.
     * @param name Name of subfield.
     */
    void _extractSubfield(const Field& field,
                          const char* name);

    /** Get scatter for given context.
     *
     * @param context Context for scatter.
     * @param createOk If true, okay to create new scatter for
     * context, if false will throw std::runtime_error if scatter for
     * context doesn't already exist.
     */
    ScatterInfo& _getScatter(const char* context,
                             const bool createOk=false);

    /** Get scatter for given context.
     *
     * @param context Context for scatter.
     */
    const ScatterInfo& _getScatter(const char* context) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    subfields_type _subfields; ///< Map of subfields in field.
    std::string _label; ///< Label for field.
    bool _dimsOkay; ///< Ok to replace nondimensionalized values with dimensionalized values.

    const Mesh& _mesh; ///< Mesh associated with section.
    scatter_map_type _scatters; ///< Collection of scatters.
    PetscDM _dm; ///< Manages the PetscSection.
    PetscVec _localVec; ///< Local PETSc vector.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Field(const Field&); ///< Not implemented
    const Field& operator=(const Field&); ///< Not implemented

}; // Field

#include "Field.icc" // inline methods

#endif // pylith_topology_field_hh

// End of file
