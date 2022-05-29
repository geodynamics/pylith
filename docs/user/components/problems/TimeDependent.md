# TimeDependent

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.TimeDependent`
:Journal name: `timedependent`

Static, quasistatic, or dynamic time-dependent problem.

Implements `Problem`.

## Pyre Facilities

* `bc`: Boundary conditions.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, bc
* `defaults`: Default options for problem.
  - **current value**: 'problem_defaults', from {default}
  - **configurable as**: problem_defaults, defaults
* `gravity_field`: Database used for gravity field.
  - **current value**: 'nullcomponent', from {default}
  - **configurable as**: nullcomponent, gravity_field
* `ic`: Initial conditions.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, ic
* `interfaces`: Interior surfaces with constraints or constitutive models.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, interfaces
* `materials`: Materials in problem.
  - **current value**: 'homogeneous', from {default}
  - **configurable as**: homogeneous, materials
* `normalizer`: Nondimensionalizer for problem.
  - **current value**: 'nondimelasticquasistatic', from {default}
  - **configurable as**: nondimelasticquasistatic, normalizer
* `petsc_defaults`: Flags controlling which default PETSc options to use.
  - **current value**: 'petscdefaults', from {default}
  - **configurable as**: petscdefaults, petsc_defaults
* `progress_monitor`: Simple progress monitor via text file.
  - **current value**: 'progressmonitortime', from {default}
  - **configurable as**: progressmonitortime, progress_monitor
* `solution`: Solution field for problem.
  - **current value**: 'solution', from {default}
  - **configurable as**: solution
* `solution_observers`: Observers (e.g., output) for solution.
  - **current value**: 'singlesolnobserver', from {default}
  - **configurable as**: singlesolnobserver, solution_observers

## Pyre Properties

* `end_time`=\<dimensional\>: End time for problem.
  - **default value**: 3.15576e+06*s
  - **current value**: 3.15576e+06*s, from {default}
  - **validator**: (greater than or equal to 0*s)
* `formulation`=\<str\>: Formulation for equations.
  - **default value**: 'quasistatic'
  - **current value**: 'quasistatic', from {default}
  - **validator**: (in ['quasistatic', 'dynamic', 'dynamic_imex'])
* `initial_dt`=\<dimensional\>: Initial time step.
  - **default value**: 3.15576e+07*s
  - **current value**: 3.15576e+07*s, from {default}
  - **validator**: (greater than 0*s)
* `max_timesteps`=\<int\>: Maximum number of time steps.
  - **default value**: 20000
  - **current value**: 20000, from {default}
  - **validator**: (greater than 0)
* `notify_observers_ic`=\<bool\>: Notify observers of solution with initial conditions.
  - **default value**: False
  - **current value**: False, from {default}
* `solver`=\<str\>: Type of solver to use ['linear', 'nonlinear'].
  - **default value**: 'linear'
  - **current value**: 'linear', from {default}
  - **validator**: (in ['linear', 'nonlinear'])
* `start_time`=\<dimensional\>: Start time for problem.
  - **default value**: 0*s
  - **current value**: 0*s, from {default}

## Example

Example of setting `TimeDependent` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Set boundary conditions, faults, and materials
bc = [boundary_xpos, boundary_xneg]
interfaces = [san_andreas, hayward]
materials = [crust, mantle]

# Create an initial condition over the domain
ic = [domain]

# Turn on gravitational body forces
gravity_field = spatialdata.spatialdb.GravityField

# Set the normalizer for nondimensionalizing the problem
normalizer = spatialdata.units.NondimElasticQuasistatic

# Set the subfields in the solution
solution = = pylith.problems.SolnDispLagrange

# Output the solution for the domain and ground surface
solution_observers = [domain, ground_surface]

# Use the quasistatic formulation, linear solver, and set appropriate default solver settings.
formulation = quasistatic
solver = linear

# Use a maximum of 20 time steps to simulation from -0.5 years to 2.0 years with an initial time step of 0.5 years.
# The first time step will compute the solution at time 0.
start_time = -0.5*year
end_time = 2.0*year
initial_dt = 0.5*year
max_timesteps = 20

[pylithapp.greensfns.petsc_defaults]
solver = True
monitors = True
:::

