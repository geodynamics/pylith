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
 * @file libsrc/topology/FieldBase.hh
 *
 * @brief Basic information related to a vector field over the
 * vertices or cells of a finite-element mesh.
 */

#if !defined(pylith_topology_fieldbase_hh)
#define pylith_topology_fieldbase_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/array.hh" // HASA string_vector

#include <string> // USES std::string
#include <map> // USES std::map

// FieldBase ------------------------------------------------------------
/** @brief Basic information related to a vector field over the
 * vertices or cells of a finite-element mesh.
 */
class pylith::topology::FieldBase {
    // PUBLIC ENUMS /////////////////////////////////////////////////////////
public:

    enum VectorFieldEnum {
        SCALAR=0, ///< Scalar.
        VECTOR=1, ///< Vector.
        TENSOR=2, ///< Tensor.
        OTHER=3, ///< Not a scalar, vector, or tensor.
        MULTI_SCALAR=4, ///< Scalar at multiple points.
        MULTI_VECTOR=5, ///< Vector at multiple points.
        MULTI_TENSOR=6, ///< Tensor at multiple points.
        MULTI_OTHER=7, ///< Not a scalar, vector, or tensor at multiple points.
    }; // VectorFieldEnum

    enum SpaceEnum {
        POLYNOMIAL_SPACE=0, ///< Polynomial finite-element space.
        POINT_SPACE=1, ///< Point finite-element space.
    }; // SpaceEnum

    enum CellBasis {
        SIMPLEX_BASIS=1, ///< Simplex basis functions.
        TENSOR_BASIS=2, ///< Tensor product basis functions.
        DEFAULT_BASIS=10, ///< Use default for cell type.
    }; // CellBasis

    // PUBLIC TYPEDEF ///////////////////////////////////////////////////////
public:

    /// Function prototype for validator functions.
    typedef const char* (*validatorfn_type)(const PylithReal);

    // PUBLIC STRUCTS ///////////////////////////////////////////////////////
public:

    struct Description {
        std::string label; ///< Label for field.
        std::string alias; ///< User-defined alias for field.
        VectorFieldEnum vectorFieldType; ///< Type of vector field.
        pylith::string_vector componentNames; ///< Names of components.
        size_t numComponents; ///< Number of components.
        PylithReal scale; ///< Dimension scale associated with values.
        validatorfn_type validator; ///< Validator for values in field;
        bool hasHistory; ///< Has subfields with history, i.e., state variables.
        size_t historySize; ///< Number of points in time history (currently only).

        Description(void) :
            label(""),
            alias(""),
            vectorFieldType(OTHER),
            numComponents(0),
            scale(1.0),
            validator(NULL),
            hasHistory(false),
            historySize(0) {}


        Description(const char* labelValue,
                    const char* aliasValue,
                    const pylith::string_vector& componentNamesValue,
                    const size_t numComponentsValue=0,
                    const VectorFieldEnum vectorFieldTypeValue=SCALAR,
                    const PylithReal scaleValue=1.0,
                    const validatorfn_type validatorValue=NULL,
                    bool isFaultOnlyValue=false,
                    bool hasHistoryValue=false,
                    const size_t historySizeValue=0) :
            label(labelValue),
            alias(aliasValue),
            vectorFieldType(vectorFieldTypeValue),
            componentNames(componentNamesValue),
            numComponents(numComponentsValue),
            scale(scaleValue),
            validator(validatorValue),
            hasHistory(hasHistoryValue),
            historySize(historySizeValue) {}


    }; // Description

    struct Discretization {
        int basisOrder; ///< Order of basis functions.
        int quadOrder; ///< Order of quadrature scheme.
        int dimension; ///< Dimension of point(s) for discretization.
        int components; ///< Number of components for field.
        bool isFaultOnly; ///< Subfield is limited to fault degrees of freedom.
        CellBasis cellBasis; ///< Cell basis (simplex, tensor, default).
        SpaceEnum feSpace; ///< Finite-element space.
        bool isBasisContinuous; ///< Is basis continuous?

        Discretization(const int basisOrderValue=1,
                       const int quadOrderValue=1,
                       const int dimensionValue=-1,
                       const int componentsValue=-1,
                       bool isFaultOnlyValue=false,
                       const CellBasis cellBasisValue=DEFAULT_BASIS,
                       const SpaceEnum feSpaceValue=POLYNOMIAL_SPACE,
                       const bool isBasisContinuousValue=true) :
            basisOrder(basisOrderValue),
            quadOrder(quadOrderValue),
            dimension(dimensionValue),
            components(componentsValue),
            isFaultOnly(isFaultOnlyValue),
            cellBasis(cellBasisValue),
            feSpace(feSpaceValue),
            isBasisContinuous(isBasisContinuousValue)
        {}


        bool operator==(const Discretization rhs) const {
            if (basisOrder         != rhs.basisOrder) {return false;}
            if (quadOrder          != rhs.quadOrder) {return false;}
            if (dimension          != rhs.dimension) {return false;}
            if (components         != rhs.components) {return false;}
            if (isFaultOnly  != rhs.isFaultOnly) {return false;}
            if (cellBasis          != rhs.cellBasis) {return false;}
            if (feSpace            != rhs.feSpace) {return false;}
            if (isBasisContinuous  != rhs.isBasisContinuous) {return false;}
            // return true;
            return false;
        } // operator=

        bool operator<(const Discretization rhs) const {
            return true;
            if (basisOrder < rhs.basisOrder) {return true;}
            if (basisOrder == rhs.basisOrder) {
                if (quadOrder < rhs.quadOrder) {return true;}
                if (quadOrder == rhs.quadOrder) {
                    if (dimension < rhs.dimension) {return true;}
                    if (dimension == rhs.dimension) {
                        if (components < rhs.components) {return true;}
                        if (components == rhs.components) {
                            if (isFaultOnly < rhs.isFaultOnly) {return true;}
                            if (isFaultOnly == rhs.isFaultOnly) {
                                if (cellBasis < rhs.cellBasis) {return true;}
                                if (cellBasis == rhs.cellBasis) {
                                    if (feSpace < rhs.feSpace) {return true;}
                                    if (feSpace == rhs.feSpace) {
                                        if (isBasisContinuous  < rhs.isBasisContinuous) {return true;}
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return false;
        } // operator<

    }; // Discretization

    /// Mapping from field name to discretization (intended for subfields).
    typedef std::map<std::string, pylith::topology::FieldBase::Discretization> discretizations_map;

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    FieldBase(void); ///< Default constructor.
    ~FieldBase(void); ///< Default destructor.

    /** Get string associated with vector field type.
     *
     * @param value Vector field type.
     * @returns String associated with vector field type.
     */
    static
    const char*vectorFieldString(VectorFieldEnum value);

    /** Get string associated with vector field type.
     *
     * @param value String associated with vector field type.
     * @returns Vector field type.
     */
    static
    VectorFieldEnum parseVectorFieldString(const char* value);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FieldBase(const FieldBase&); ///< Not implemented
    const FieldBase& operator=(const FieldBase&); ///< Not implemented

}; // FieldBase

#endif // pylith_topology_fieldbase_hh

// End of file
