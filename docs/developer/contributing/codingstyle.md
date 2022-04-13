# Coding Style

There are a number of standard coding styles for programming languages, notably PEP8 for Python. For PyLith, we try to be consistent in naming conventions across Python and C++ while following a subset of the conventions used in PETSc and PEP8 with documentation styles consistent with Doxygen.

:::{important}
We use 4 spaces for indentation. Configure your editor to use spaces instead of tabs.
:::

## General guidelines

* Naming conventions
  * Use self-documenting names.
  * Avoid single letter variables. Choose names that are readily found via searches across single or multiple files (e.g., grep).
  * Class names are generally nouns and methods are verbs.
  * Class names use upper camel case, e.g., `TimeDependent`.
  * Public method names use camel case, e.g., `computeRHSResidual()`.
  * Protected and private method names use camel case preceded by an underscore, e.g., `_setFEKernelsRHSResidual()`.
  * In C++ data members are private and use camel case preceded by an underscore, e.g., `_gravityField`.
  * In Python data members are public and use camel case, e.g., `self.gravityField`.
  * Local variables use camel case, e.g., `numIntegrators`.
* Comments
  * List authors, copyright, and license info at the very beginning of every file.
  * For every class method, describe its function and include a description for each argument. For Python this is done in the docstring of the method, and for C++ this is done in a doxygen style comment immediately before the method declaration in the header file.
  * Document nontrivial algorithms and any assumptions.
* Error checking
  * PyLith should never crash without an error message.
  * All user errors should be trapped as early as possible and reported with an informative error message by throwing an appropriate exception. If possible, suggest ways to correct the error.
  * Messages for internal errors should indicate the location in the code where the error was trapped.
  * All pointers should be checked for `NULL` values before use. Usually we use `assert()` to do this.
  * Check the return values for all calls to functions in external libraries. For PETSc functions, we use `PYLITH_CHECK_ERROR(returnValue)`.
* Testing
  * All C++ methods should be covered by unit tests.
  * All governing equations should be covered by Method of Manufactured solution tests.
  * All functionality should be covered by full-scale tests.

## Formatting source files

We use `autopep8` and `uncrustify` to format Python and C//C++ source files, respectively. The corresponding configuration files are `autopep8.cfg` and `uncrustify.cfg` in the `developer` directory. The Python script `developer/format_source.py` is a handy utility for calling `autopep8` and `uncrustify` with the appropriate arguments and formatting multiple files.

:::{tip}
We highly recommend using an integrated development environment with `uncrustify` and `autopep8` to automatically format all C/C++ and Python source files.
:::

```{code-block} bash
---
caption: Formatting Python and C++ source code using `format\_source.py`. `autopep8` and `uncrustify` must be in the current path.
---
# Format a C++ file.
developer/format_source.py --cplusplus=libsrc/pylith/materials/Material.cc

# Format all C++ files in the 'libsrc/pylith/materials' directory.
developer/format_source.py --cplusplus=libsrc/pylith/materials/*.cc

# Format a Python file.
developer/format_source.py --python=pylith/materials/Material.py

# Format all Python files in the 'pylith/materials' directory.
developer/format_source.py --python=pylith/materials/*.py
```


## Error Checking

Our philosophy is that PyLith should never crash without an error message. If it encounters a fatal error, then it should generate an appropriate error message and abort. In C++ we throw `std::runtime_error` exceptions for errors resulting from user input and `std::logic_error` exceptions for internal inconsistencies or logic errors. In Python we use standard exception objects.

Additional protections against crashing include: using asserts to verify pointers are non-null before using them and using the `PYLITH_CHECK_ERROR` macro to check the return value after *every* call to a PETSc function.

```{code-block} c++
---
caption: Example of using `assert()`
emphasize-lines: 1
---

assert(_solution); // Verify _solution is not NULL.

// Initialize integrators.
const size_t numIntegrators = _integrators.size();
for (size_t i = 0; i < numIntegrators; ++i) {
    assert(_integrators[i]); // Verify _integrators[i] is not NULL.
    _integrators[i]->initialize(*_solution);
} // for  
```

:::{tip}
When we build the code for production runs, we usually configure with `CPPFLAGS=-DNDEBUG` to remove `assert()` calls.
:::

```{code-block} c++
---
caption: Example of using `PYLITH_CHECK_ERROR` macro.
---
PetscErrorCode err = TSGetTimeStep(ts, &dt); PYLITH_CHECK_ERROR(err);
```

In combination with the above procedures, we also make use of the Pyre journals to display warnings and errors to facilitate debugging. The journals provide the file name and line number along with the message. By default, Pyre journals for errors are turned on and those for warnings and debugging are turned off.

```{code-block} c++
---
caption: Example of using Pyre journals and standard exceptions.
emphasize-lines: 9, 12-13
---
switch (bitUse) {
case 0x1:
    _bcKernel = pylith::fekernels::TimeDependentFn::initial_scalar;
    break;
case 0x2:
    _bcKernel = pylith::fekernels::TimeDependentFn::rate_scalar;
    break;
case 0x0:
    PYLITH_COMPONENT_WARNING("Dirichlet BC provides no constraints.");
    break;
default:
    PYLITH_COMPONENT_LOGICERROR("Unknown combination of flags for Dirichlet BC terms "
        << "(useInitial="<<_useInitial<<", useRate="<<_useRate<<").");
} // switch
```

## C/C++ style

### Object Declaration Files

Object declaration (header) files use the `.hh` suffix. C header files use the `.h` suffix. The following code excerpt demonstrates the conventions we use in formatting header files and including comments.

:::{important}
*All* declarations of class methods should include a description of what the method does and a description of each argument and the return value if it is not void.
:::

```{code-block} c++
---
caption: Sample C++ declaration (header) file.
---
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, The State University of New York at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2018 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/* ANNOTATION: Next list the filename with the relative path of the file
 * along with a brief description.
 */  

/**
 * @file libsrc/problems/Problem.hh
 *
 * @brief C++ object that manages formulating the equations.
 */

// ANNOTATION: Use full namespace in header guard.
#if !defined(pylith_problems_problem_hh)
#define pylith_problems_problem_hh

/* ANNOTATION: Includes.
 *
 * 1. Header file for forward declaration of class
 * 2. Header file for parent classes
 * 3. Local forward declarations
 * 4. Local header files
 * 5. Other header files
 * 6. System header files
 */
#include "problemsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/feassemble/feassemblefwd.hh" // HASA IntegratorPointwise
#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field
#include "pylith/meshio/meshiofwd.hh" // HASA OutputManager
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA GravityField

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat

#include "pylith/utils/array.hh" // HASA std::vector


/* ANNOTATION: Provide description of class.
 */

/** Reform the Jacobian and residual for the problem.
 *
 * We cast the problem in terms of F(t,s,\dot{s}) = G(t,s), s(t0) = s0.
 *
 * In PETSc time stepping (TS) notation, G is the RHS, and F is the I
 * function (which we call the LHS).
 *
 */
class pylith::problems::Problem : public pylith::utils::PyreComponent {
  /* ANNOTATION: Order of declarations is:
  *
  * 1. Friend classes
  * 2. Public enums
  * 3. Public structs
  * 4. Public methods
  * 5. Protected methods
  * 6. Private methods
  * 7. Protected members
  * 8. Private members
  * 9. Methods not implemented
  *
  * Within each group, we generally order the methods by:
  * 1. Constructor/destructors
  * 2. Accessors
  * 3. Other methods
  *
  * Use the full namespace when declaring data members and method
  * arguments to avoid ambiguity.
  *
  * Method arguments are listed one per line.
  *
  * Include opening braces at the end of a line. Use a comment to
  * document all closing braces.
  *
  * Before every member method, describe what the method does and
  * include a description for every argument. We use Doxygen
  * syntax.
  */

  friend class TestProblem;   // unit testing

    // PUBLIC ENUM //////////////////////////////////////////////////////////
public:

    enum SolverTypeEnum {
        LINEAR, // Linear solver.
        NONLINEAR, // Nonlinear solver.
    };   // SolverType


    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // Constructor
    Problem(void);

    /// Destructor
    virtual ~Problem(void);

    /* We call the deallocate method before calling PetscFinalize() to
     * deallocate any memory allocated using PETSc. In general, the
     * destructor will simply call deallocate().
     */

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set solver type.
     *
     * @param[in] value Solver type.
     */
    void solverType(const SolverTypeEnum value);

    /** Get solver type.
     *
     * @returns Solver type.
     */
    SolverTypeEnum solverType(void) const;

    /** Set manager of scales used to nondimensionalize problem.
     *
     * @param[in] dim Nondimensionalizer.
     */
    void normalizer(const spatialdata::units::Nondimensional& dim);

    /** Set gravity field.
     *
     * @param[in] g Gravity field.
     */
    void gravityField(spatialdata::spatialdb::GravityField* const g);

    /** Set solution field.
     *
     * @param[in] field Solution field.
     */
    void solution(pylith::topology::Field* field);

    /** Set handles to integrators.
     *
     * @param[in] integratorArray Array of integrators.
     * @param[in] numIntegrators Number of integrators.
     */
    void integrators(pylith::feassemble::IntegratorPointwise* integratorArray[],
                     const int numIntegrators);

    /** Do minimal initialization.
     *
     * @param mesh Finite-element mesh.
     */
    virtual
    void preinitialize(const pylith::topology::Mesh& mesh);

    /** Verify configuration.
     *
     * @param[in] materialIds Array of material ids.
     * @param[in] numMaterials Size of array (number of materials).
     *
     */
    virtual
    void verifyConfiguration(int* const materialIds,
                             const int numMaterials) const;

    /** Set solution values according to constraints (Dirichlet BC).
     *
     * @param[in] t Current time.
     * @param[in] solutionVec PETSc Vec with current global view of solution.
     * @param[in] solutionDotVec PETSc Vec with current global view of time derivative of solution.
     */
    void setSolutionLocal(const PylithReal t,
                          PetscVec solutionVec,
                          PetscVec solutionDotVec);

    /** Compute RHS residual, G(t,s) and assemble into global vector.
     *
     * @param[out] residualVec PETSc Vec for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     */
    void computeRHSResidual(PetscVec residualVec,
                            const PetscReal t,
                            const PetscReal dt,
                            PetscVec solutionVec);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    /* ANNOTATION: Use pointers to hide implementation details and speed up compilation.
     */

    pylith::topology::Field* _solution;   ///< Handle to solution field.
    pylith::topology::Field* _residual; ///< Handle to residual field.

    spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalization of scales.
    spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.
    std::vector<pylith::feassemble::IntegratorPointwise*> _integrators;   ///< Array of integrators.
    SolverTypeEnum _solverType;   ///< Problem (solver) type.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    /* ANNOTATION: Declare expensive, fragile copy methods private, so using them
     * fails at compile time with an error.
     */

    Problem(const Problem&);   ///< Not implemented
    const Problem& operator=(const Problem&);   ///< Not implemented

}; // Problem

#endif // pylith_problems_problem_hh


// End of file
```

### Object Implementation Files

Object implementation files use the `.cc` suffix. Inline implementation files use the `.icc` suffix and are included from the definition files. C implementation files use the `.c` suffix.

To facilitate debugging and error messages, we use the following
macros:

`PYLITH_METHOD_BEGIN`
: This macro allows line numbers of source files to be included in PETSc error messages. Use this macro at the beginning of all methods using any PETSc routines as well as most other methods. We don't use this macro in destructors because many of them are called *after* `PetscFinalize`. We also do not use this macro in trivial or inline methods that do not call any PETSc routines.

`PYLITH_METHOD_END`
: Use the macro at the end of all methods that begin with `PYLITH_METHOD_BEGIN` and return void.

`PYLITH_RETURN_END`
: Use this macro at the end of all methods that begin with `PYLITH_METHOD_BEGIN` and return non-void values.

`PYLITH_CHECK_ERROR`
: Use this macro after *every* call to a PETSc function to check the return value.

`PYLITH_JOURNAL_DEBUG`
: Use this macro immediately after `PYLITH_METHOD_BEGIN` in methods of all objects that inherit from `GenericComponent`.

`PYLITH_COMPONENT_DEBUG`
: Use this macro immediately after `PYLITH_METHOD_BEGIN` in methods of all objects that inherit from `PyreComponent`.  Non-abstract classes should call `PyreComponent::setName()` in the constructor. We recommend using a static data member for the name with the lowercase name matching the Pyre component, e.g., "timedependent" for the C++ `TimeDependent` object.

```{code-block} c++
---
caption: Sample C++ definition (implementation) file.
---
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, The State University of New York at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/* ANNOTATION: Order of including header files is:
 * 1. portinfo: (stuff from configure)
 * 2. Header file for this class.
 * 3. Header files for local classes.
 * 4. Header files for classes in other libraries.
 * 5. Standard header files.
 *
 * List why each header file is included (USES/HASA/HOLDSA).
 *
 * Order of method implementations should match header file.
 *
 * For each method:
 *   1. Return values should go on previous line.
 *   2. Put each method argument on a separate line.
 */
#include <portinfo>

#include "Problem.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/feassemble/IntegratorPointwise.hh" // USES IntegratorPointwise

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem() :
    _solution(NULL),
    _solutionDot(NULL),
    _residual(NULL),
    _jacobianLHSLumpedInv(NULL),
    _normalizer(NULL),
    _gravityField(NULL),
    _integrators(0),
    _constraints(0),
    _outputs(0),
    _solverType(LINEAR) {}

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Problem::~Problem(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Problem::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _solution = NULL; // Held by Python. :KLUDGE: :TODO: Use shared pointer.
    delete _solutionDot; _solutionDot = NULL;
    delete _residual; _residual = NULL;
    delete _jacobianLHSLumpedInv; _jacobianLHSLumpedInv = NULL;
    delete _normalizer; _normalizer = NULL;
    _gravityField = NULL; // Held by Python. :KLUDGE: :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set problem type.
void
pylith::problems::Problem::solverType(const SolverTypeEnum value) {
    PYLITH_COMPONENT_DEBUG("Problem::solverType(value="<<value<<")");

    _solverType = value;
} // solverType

// ----------------------------------------------------------------------
// Get problem type.
pylith::problems::Problem::SolverTypeEnum
pylith::problems::Problem::solverType(void) const {
    return _solverType;
} // solverType

// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Problem::normalizer(const spatialdata::units::Nondimensional& dim) {
    PYLITH_COMPONENT_DEBUG("Problem::normalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer

// ----------------------------------------------------------------------
// Set gravity field.
void
pylith::problems::Problem::gravityField(spatialdata::spatialdb::GravityField* const g) {
    PYLITH_COMPONENT_DEBUG("Problem::gravityField(g="<<typeid(*g).name()<<")");

    _gravityField = g;
} // gravityField

// ----------------------------------------------------------------------
// Set solution field.
void
pylith::problems::Problem::solution(pylith::topology::Field* field)
{ // solution
    PYLITH_COMPONENT_DEBUG("Problem::solution(field="<<typeid(*field).name()<<")");

    _solution = field;
} // solution

// ----------------------------------------------------------------------
// Set integrators over the mesh.
void
pylith::problems::Problem::integrators(pylith::feassemble::IntegratorPointwise* integratorArray[],
                                       const int numIntegrators) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::integrators("<<integratorArray<<", numIntegrators="<<numIntegrators<<")");

    assert( (!integratorArray && 0 == numIntegrators) || (integratorArray && 0 < numIntegrators) );

    _integrators.resize(numIntegrators);
    /* Declare loop variables inline. Always use braces at begin/end
     * of if and for statements.
     */
    for (int i = 0; i < numIntegrators; ++i) {
        _integrators[i] = integratorArray[i];
    } // for

    PYLITH_METHOD_END;
} // integrators

// ----------------------------------------------------------------------
// Do minimal initialization.
void
pylith::problems::Problem::preinitialize(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::preinitialzie(mesh="<<typeid(mesh).name()<<")");

    assert(_normalizer);

    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->normalizer(*_normalizer);
        _integrators[i]->gravityField(_gravityField);
    } // for

    PYLITH_METHOD_END;
} // preinitialize

// ----------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::Problem::verifyConfiguration(int* const materialIds,
                                               const int numMaterials) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(materialIds="<<materialIds<<", numMaterials="<<numMaterials<<")");

    assert(_solution);

    // Check to make sure material-id for each cell matches the id of a material.
    pylith::topology::MeshOps::checkMaterialIds(_solution->mesh(), materialIds, numMaterials);

    // Check to make sure integrators are compatible with the solution.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->verifyConfiguration(*_solution);
    } // for

    PYLITH_METHOD_END;
}  // verifyConfiguration

// ----------------------------------------------------------------------
// Set solution values according to constraints (Dirichlet BC).
void
pylith::problems::Problem::setSolutionLocal(const PylithReal t,
                                            PetscVec solutionVec,
                                            PetscVec solutionDotVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolutionLocal(t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<")");

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterVectorToLocal(solutionVec);

    if (solutionDotVec) {
        if (!_solutionDot) {
            _solutionDot = new pylith::topology::Field(_solution->mesh());
            _solutionDot->cloneSection(*_solution);
            _solutionDot->setLabel("solutionDot");
        } // if
        _solutionDot->scatterVectorToLocal(solutionDotVec);
    } // if

    PYLITH_METHOD_END;
} // setSolutionLocal

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::problems::Problem::computeRHSResidual(PetscVec residualVec,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeRHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(_solution);

    // Update PyLith view of the solution.
    PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum residual contributions across integrators.
    _residual->zeroLocal();
    const size_t numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeRHSResidual(_residual, t, dt, *_solution);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0); PYLITH_CHECK_ERROR(err); // Move to TSComputeIFunction()?
    _residual->scatterLocalToVector(residualVec, ADD_VALUES);

    PYLITH_METHOD_END;
} // computeRHSResidual

// End of file
```

## Python style

```{code-block} python
---
caption: Sample Python source code.
---
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, The State University of New York at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# ANNOTATION: Next list the filename with the relative path of the file along with a
# brief description.
#  
# @file pylith/problems/Problem.py
#
# @brief Python abstract base class for crustal dynamics problems.

# ANNOTATION: Order of imports should be
# 1. standard modules
# 2. other modules
# 3. PyLith modules (parent classes first).
#
# The SWIG interface to the C++ object is imported as ModuleProblem in
# order to avoid clashing with the Python Problem class.

from pylith.utils.PetscComponent import PetscComponent
from .problems import Problem as ModuleProblem

# ITEM FACTORIES ///////////////////////////////////////////////////////

# ANNOTATION: Define any factory methods needed in the Pyre inventory.
def materialFactory(name):
    """
    Factory for material items.
    """
    from pyre.inventory import facility
    from pylith.materials.ElasticityPlaneStrain import ElasticityPlaneStrain
    return facility(name, family="material", factory=IsotropicLinearElasticityPlaneStrain)

class Problem(PetscComponent, ModuleProblem):
    """Python abstract base class for crustal dynamics problems.

    FACTORY: problem. List the factory if one exists.
    """

    import pyre.inventory

    # ANNOTATION: Usually, we put all arguments on a single line. If the
    # function call is really long, break the arguments into
    # logical pieces, usually one argument per line.
    solverTypeStr = pyre.inventory.str(
        "solver",
        default="linear",
        validator=pyre.inventory.choice(["linear", "nonlinear"])
    )
    solverTypeStr.meta['tip'] = "Type of solver to use ['linear', 'nonlinear']."

    from Solution import Solution
    solution = pyre.inventory.facility("solution", family="solution", factory=Solution)
    solution.meta['tip'] = "Solution field for problem."

    from pylith.materials.Homogeneous import Homogeneous
    materials = pyre.inventory.facilityArray(
        "materials",
        itemFactory=materialFactory,
        factory=Homogeneous
    )
    materials.meta['tip'] = "Materials in problem."


    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="problem"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="problem")

        # ANNOTATION: Initialize all data members not in the inventory.
        self.mesh = None
        return

    def preinitialize(self, mesh):
        """Do minimal initialization.
        """
        # ANNOTATION: On process 0 only, print progress information to info journal.
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization before verifying configuration.")

        # ANNOTATION: Pass information to corresponding C++ object.
        #
        # For all calls to the C++ interaface, call using the
        # class name and pass self as an argument to make it clear
        # that this is calling a C++ method and not a Python method.
        ModuleProblem.setIdentifier(self, self.aliases[-1])

        if self.solverTypeStr == "linear":
            ModuleProblem.solverType(self, ModuleProblem.LINEAR)
        elif self.solverTypeStr == "nonlinear":
            ModuleProblem.solverType(self, ModuleProblem.NONLINEAR)
        else:
            raise ValueError("Unknown solver type '%s'." % self.solverTypeStr)

        # Do minimal setup of solution.
        self.solution.preinitialize(mesh, self.normalizer)
        ModuleProblem.solution(self, self.solution.field)

        # Preinitialize materials
        for material in self.materials.components():
            material.preinitialize(mesh)

        ModuleProblem.preinitialize(self, mesh)
        return

    def initialize(self):
        """Initialize integrators and constraints.
        """
        # On process 0 only, print progress information to info journal.
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Initializing problem.")

        ModuleProblem.initialize(self)
        return

    def run(self, app):
        """Solve the problem.
        """
        # Generate error if method is not implemented in child class.
        raise NotImplementedError("run() not implemented.")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set data members based using inventory.
        """
        PetscComponent._configure(self)

        return

# FACTORIES ////////////////////////////////////////////////////////////

# ANNOTATION: Define any facility factories.
def problem():
    """
    Factory associated with Problem.
    """
    return Problem()


# End of file
```

% End of file
