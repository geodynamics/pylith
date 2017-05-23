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
 * 1. Set the number of subfields and the number of components per subfield.
 * 2. Set the chart.
 * 3. Set the DOF (fiberdim) for each point in the chart.
 * 4. Set any constraints (known DOF).
 * 5. Allocate.
 * 6. Set the indices for the constrained DOF.
 *
 * For local fields associated with parameters, we do not have steps 4
 * and 6. Step 2-4 and 6 are automatically handled by the PETSc DS if
 * DMAddBoundary is called.
 */
class pylith::topology::Field : public FieldBase, public pylith::utils::GenericComponent {
    friend class FieldQuery;   // Fill field using data.

    friend class TestFieldMesh;   // unit testing
    friend class TestFieldSubMesh;   // unit testing

    // PUBLIC STRUCTS ///////////////////////////////////////////////////////
public:

    /// Subfield auxiliary information.
    struct SubfieldInfo {
        Description description; ///< Description for subfield.
        Discretization fe; ///< Discretization information for subfield.
        int index; ///< Index of subfield in field.
        PetscDM dm; ///< PETSc DM associated with subfield.
    };   // SubfieldInfo

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param mesh Finite-element mesh.
     */
    Field(const Mesh& mesh);

    #if 0
    /** Constructor with mesh, DM, and metadata
     *
     * @param mesh Finite-element mesh.
     * @param dm PETSc dm for field.
     * @param label Label of field.
     */
    Field(const Mesh& mesh,
          PetscDM dm,
          const char* label);
          #endif

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

    /** Has section been setup?
     *
     * @returns True if section has been setup.
     */
    bool hasSection(void) const; // :TODO: @brad Remove, obsolete.

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

    /// Set chart for solution.
    void setupSolnChart(void); // :TODO: @brad Remove, obsolete.

    /** Set default DOF for solution.
     *
     * @param fiberDim Total number of components in solution.
     * @param subfieldName Name of subfield for DOF.
     */
    void setupSolnDof(const int fiberDim, // :TODO: @brad Remove obsolte.
                      const char* subfieldName="displacement");

    /** Create PETSc section and set chart and fiber dimesion for a list
     * of points.
     *
     * @param pStart First point
     * @param pEnd Upper bound for points
     * @param dim Fiber dimension for section.
     *
     * @note Don't forget to call label(), especially if reusing a field.
     */
    void newSection(const PetscInt pStart, // :TODO: @brad Remove obsolete.
                    const PetscInt pEnd,
                    const int fiberDim);

    /** Create PETSc section and set chart and fiber dimesion for a list
     * of points.
     *
     * @param points Points over which to define section.
     * @param dim Fiber dimension for section.
     *
     * @note Don't forget to call label(), especially if reusing a field.
     */
    void newSection(const int_array& points, // :TODO: @brad Remove obsolete.
                    const int fiberDim);

    /** Create PETSc section and set chart and fiber dimesion for a list
     * of points.
     *
     * @param points Points over which to define section.
     * @param num The number of points
     * @param dim Fiber dimension for section.
     *
     * @note Don't forget to call label(), especially if reusing a field.
     */
    void newSection(const PetscInt *points,
                    const PetscInt num,
                    const int fiberDim);

    /** Create PETSc section and set chart and fiber dimesion.
     *
     * @param domain Type of points over which to define section.
     * @param dim Fiber dimension for section.
     * @param stratum Stratum depth (for vertices) and height (for cells).
     *
     * @note Don't forget to call label(), especially if reusing a field.
     */
    void newSection(const DomainEnum domain, // :TODO: @brad Remove obsolete.
                    const int fiberDim,
                    const int stratum=0);

    /** Create section using src field as template with given fiber dimension.
     *
     * @param sec Field defining layout.
     * @param fiberDim Fiber dimension.
     *
     * @note Don't forget to call label(), especially if reusing a field.
     */
    void newSection(const Field& src,
                    const int fiberDim);

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
     * @param[in] name Name of subfield.
     * @param[in] fieldType Type of vector field.
     * @param[in] components Array of names of field components.
     * @param[in] numComponents Size of array.
     * @param[in] scale Dimensional scale associated with field.
     * @param[in] basisOrder Order of basis functions for discretization.
     * @param[in] quadOrder Order of numerical quadrature for discretization.
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space (POLYNOMIAL_SPACE or POINT_SPACE).
     */
    void subfieldAdd(const char *name,
                     const VectorFieldEnum fieldType,
                     const char* components[],
                     const int numComponents,
                     const double scale,
                     const int basisOrder,
                     const int quadOrder,
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

    /** Setup sections for subfields.
     *
     * Should be preceded by calls to subfieldAdd() and followed by calls to subfieldSetDof().
     */
    void subfieldsSetup(void);

    /** Convenience method for setting number of DOF (fiberdim) for subfield at points.
     *
     * Should be preceded by calls to subfieldAdd() and subfieldsSetup().
     *
     * @param name Name of subfield.
     * @param domain Point classification for subfield.
     * @param fiberDim Number of subfield components per point.
     */
    void subfieldSetDof(const char *name,  // :TODO: @brad Remove obsolete.
                        const pylith::topology::FieldBase::DomainEnum domain,
                        const int fiberDim);

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
     * @param label Label for output.
     */
    void view(const char* label) const;

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
     * @param mesh Mesh associated with scatter.
     * @param context Label for context associated with vector.
     */
    void createScatterWithBC(const Mesh& mesh,
                             const char* context);


    /** Create PETSc vector scatter for field. This is used to transfer
     * information from the "global" PETSc vector view to the "local"
     * PETSc section view. The PETSc vector includes constrained
     * DOF. Use createScatter() if constrained DOF should be omitted
     * from the PETSc vector.
     *
     * @param mesh Mesh associated with scatter.
     * @param labelName The name of the label defining the point set, or PETSC_NULL
     * @param labelValue The label stratum defining the point set
     * @param context Label for context associated with vector.
     */
    void createScatterWithBC(const Mesh& mesh,
                             const std::string& labelName,
                             PetscInt labelValue,
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
    };   // ScatterInfo

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

    subfields_type _subfields;   ///< Map of subfields bundled together.
    std::string _label; ///< Label for field.
    bool _dimsOkay; ///< Ok to replace nondimensionalized values with dimensionalized values.

    const Mesh& _mesh;   ///< Mesh associated with section.
    scatter_map_type _scatters;   ///< Collection of scatters.
    PetscDM _dm;   ///< Manages the PetscSection.
    PetscVec _localVec;   ///< Local PETSc vector.

    static const char* _genericComponent; ///< Name of generic component.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Field(const Field&);   ///< Not implemented
    const Field& operator=(const Field&);   ///< Not implemented

}; // Field

#include "Field.icc" // inline methods

#endif // pylith_topology_field_hh


// End of file
