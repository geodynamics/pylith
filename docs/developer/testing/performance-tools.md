(developer-performance-tools)=
# Performance Tools

:::{note}
New in v4.1.0
:::

## Overview

We use PETSc event logging to instrument the PyLith code for evaluating runtime performance.
This allows us to simultaneously evaluate the performance of the PyLith C++ code as well as the PETSc code.
In PyLith v2, we logged the floating point operations for the residual and Jacobian evaluation.
We have not yet included these in PyLith v3 and later because we now use kernels that are called many, many times and event logging in this methods would introduce overhead and reduce performance.

## Implementation of Event Logging

PyLith includes `pylith::utils::EventLogger` as a high-level interface for event logging with PETSc.
We create local `Events` classes within the C++ implementation files (`*.cc`) to create variables that hold the event ids.
This streamlines logging events within the C++ code.
We identify events using variables rather than strings, which means that typos in names are detected at compile time rather than at runtime.
Refer to the `_Problem::Events` class near the top of `libsrc/pylith/problems/Problem.cc` for an example of how we implement event logging.

We name events in PyLith C++ code following the template `PL:CLASS:METHOD` (C++ namespace style but reducing `::` to `:` to reduce the length of the string) in which we are logging method `METHOD` in C++ class `CLASS`.
In the PyLith Python code, we follow the template `PL.CLASS.METHOD`.

## Evaluating Performance

:::{tip}
We can get a general ideal of runtime performance bottlenecks with a debugging build of the code, but fine-tuning performance should be done using an optimized build.
:::

### Quick View

Run a PyLith simulation adding the command line argument `--petsc.log_view=ascii:FILENAME`.
This will dump performance information in ASCII format to `FILENAME`.

### Graphical View

Run a PyLith simulation adding the command line argument `--petsc.log_view=:FILENAME_LOGVIEW:ascii_flamegraph`. This will dump performance information to `FILENAME_LOGVIEW` that can be viewed using flame graph tools, such as [FlameGraph](https://github.com/brendangregg/FlameGraph).

#### Viewing Performance with FlameGraph

1. Run `PATH_TO_FLAMEGRAPH/flamegraph.pl FILENAME_LOGVIEW > FILENAME_PLOT.svg`
2. Load `FILENAME_PLOT.svg` into your favorite `svg` viewer, such as a web browser.

:::{note}
You can install [FlameGraph](https://github.com/brendangregg/FlameGraph) by cloning the repository or downloading and unpacking a release.
:::
