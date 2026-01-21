#!/usr/bin/env nemesis

import importlib

from matplotlib import pyplot
import numpy

import matplotlib_extras

N_UNKNOWNS = numpy.array([3878, 15128, 59696, 237104, 945008, 3773168])

# grep LocalTimes STEP.py
# grep "PL:TimeDependent:solve" STEP.py
STATS = {
    # "step01_fieldsplit_full": (6, 6, 8,),
    "step02_fieldsplit_selfp": {
        "n_iterations": (45, 104, numpy.nan),
        "total_time": (11.3959, 25.3932, numpy.nan),
        "solve_time": (4.537, 13.5346, numpy.nan),
    },
    "step03_vpbjacobi": {
        "n_iterations": (29, 42, 59, 136, 490),
        "total_time": (12.4386, 24.8782, 80.3902, 319.493, 1606.67),
        "solve_time": (4.81182, 12.4584, 46.0547, 197.797, 1138.08),
    },
    "step04_vpbjacobi_tunefine": {
        "n_iterations": (15, 19, 21, 24, 32, 56),
        "total_time": (12.0174, 24.4725, 78.0262, 296.251, 1184.62, 5087.16),
        "solve_time": (4.64184, 12.2738, 44.3372, 175.193, 715.98, 3155.8),
    },
    "step05_vpbjacobi_tunecoarse": {
        "n_iterations": (16, 20, 26, 34, 56, 58),
        "total_time": (11.9595, 25.8678, 83.454, 320.817, 1337.99, 5552.66),
        "solve_time": (4.80032, 13.6441, 50.1707, 202.285, 870.153, 3637.25),
    },
    "step06_vpbjacobi_dispscale": {
        "n_iterations": (16, 20, 25, 32, 40, 69),
        "total_time": (11.6945, 25.7422, 84.2219, 323.161, 1348.64, 5553.58),
        "solve_time": (4.60636, 13.747, 50.9383, 203.695, 837.982, 3664.83),
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
