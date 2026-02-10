#!/usr/bin/env nemesis

import importlib

from matplotlib import pyplot
import numpy

import matplotlib_extras

N_UNKNOWNS = numpy.array([3878, 15128, 59696, 237104, 945008, 3773168])

# grep LocalTimes STEP.py
# grep "PL:TimeDependent:solve" STEP.py
STATS = {
    "step02_fieldsplit_selfp": {
        "n_iterations": (1,),
        "total_time": (0,),
        "solve_time": (0,),
    },
    "step03_vpbjacobi": {
        "n_iterations": (48, numpy.nan),
        "total_time": (69.9805, numpy.nan),
        "solve_time": (33.8381, numpy.nan),
    },
    "step04_vpbjacobi_tunefine": {
        "n_iterations": (21, ),
        "total_time": (69.7465,),
        "solve_time": (34.2697,),
    },
    "step05_vpbjacobi_tunecoarse": {
        "n_iterations": (),
        "total_time": (),
        "solve_time": (),
    },
    "step06_vpbjacobi_dispscale": {
        "n_iterations": (),
        "total_time": (),
        "solve_time": (),
    },
}

LABELS = (
    "fieldsplit (selfp)",
    "vpbjacobi (defaults)",
    "vpbjacobi (tunefine)",
    "vpbjacobi (tunecoarse)",
    "vpbjacobi (dispscale)",
)

pyplot.style.use("matplotlib_extras.color-darkbg")
figure, axes = pyplot.subplots(
    ncols=2, nrows=1, figsize=(8.5, 4.0), layout="constrained"
)

ax_iters = axes[0]
ax_iters.spines[["top", "right"]].set_visible(False)
ax_iters.set_xlabel("Number of unknowns")
ax_iters.set_ylabel("Number of solver iterations")
ax_iters.set_xlim(0.8 * numpy.min(N_UNKNOWNS), 1.2 * numpy.max(N_UNKNOWNS))
ax_iters.set_ylim(10, 20000)
ax_iters.set_xscale("log")
ax_iters.set_yscale("log")

ax_time = axes[1]
ax_time.spines[["top", "right"]].set_visible(False)
ax_time.set_xlabel("Number of unknowns")
ax_time.set_ylabel("Time, s")
ax_time.set_xlim(0.8 * numpy.min(N_UNKNOWNS), 1.2 * numpy.max(N_UNKNOWNS))
ax_time.set_ylim(1, 50000)
ax_time.set_xscale("log")
ax_time.set_yscale("log")

for i_solver, sim_values in enumerate(STATS.values()):
    n_sims = len(sim_values["n_iterations"])

    line = ax_iters.plot(
        N_UNKNOWNS[:n_sims],
        sim_values["n_iterations"],
        "s-",
        lw=2,
        label=LABELS[i_solver],
    )
    if numpy.isnan(sim_values["n_iterations"])[-1]:
        i_nan = n_sims - 1
        max_y = float(ax_iters.get_ylim()[1])
        ax_iters.plot(
            N_UNKNOWNS[i_nan - 1 : i_nan + 1],
            [sim_values["n_iterations"][i_nan - 1], max_y],
            "--",
            lw=2,
            color=line[0].get_color(),
        )

    line = ax_time.plot(
        N_UNKNOWNS[:n_sims],
        sim_values["total_time"],
        "s-",
        lw=2,
        label=LABELS[i_solver] + " total",
    )
    ax_time.plot(
        N_UNKNOWNS[:n_sims],
        sim_values["solve_time"],
        "^-",
        lw=1,
        color=line[0].get_color(),
        label=LABELS[i_solver] + " solve",
    )

    ax_time.legend(loc="upper left", fontsize="small")
    figure.savefig(f"solver-performance-{i_solver}.pdf")
