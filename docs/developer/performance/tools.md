(developer-performance-tools)=
# Performance Tools

*New in v4.1.0.*

## Overview

We use PETSc event logging to instrument the PyLith code for evaluating runtime performance.
This allows us to simultaneously evaluate the performance of the PyLith C++ code as well as the PETSc code.
In PyLith v2, we logged the floating point operations for the residual and Jacobian evaluation.
We have not yet included these in PyLith v3 and later because we now use kernels that are called many, many times and event logging in this methods would introduce overhead and reduce performance.

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
