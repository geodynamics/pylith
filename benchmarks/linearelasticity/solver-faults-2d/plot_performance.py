#!/usr/bin/env nemesis

from matplotlib import pyplot
import numpy

N_UNKNOWNS = numpy.array([3878, 15128, 59696, 237104, 945008, 3773168])

# grep LocalTimes STEP.py
# grep "PL:TimeDependent:solve" STEP.py
STATS = {
    "step02_fieldsplit_selfp": {
        "n_iterations": (37, 84, 1017, numpy.nan),
        "total_time": (5.78073, 7.67613, 55.4353, numpy.nan),
        "solve_time": (1.98163, 3.51442, 48.6213, numpy.nan),
    },
    "step03_vpbjacobi": {
        "n_iterations": (25, 37, 47, 141, 157, 1182),
        "total_time": (5.8148, 7.38562, 14.615, 56.6935, 226.305, 3048.99),
        "solve_time": (2.02259, 3.06247, 8.24445, 40.2364, 171.427, 2824.93),
    },
    "step04_vpbjacobi_tunefine": {
        "n_iterations": (12, 16, 17, 20, 25, 31),
        "total_time": (5.62875, 6.98172, 13.5789, 43.8165, 165.896, 712.243),
        "solve_time": (1.68937, 2.88792, 7.32962, 27.8018, 112.093, 495.028),
    },
    # Not as good as step04
    # "step05_vpbjacobi_tunecoarse": {
    #    "n_iterations": (16, 17, 22, 28, 34, 41),
    #    "total_time": (5.5662, 7.12535, 15.5069, 50.3177, 197.339, 860.028),
    #    "solve_time": (1.80709, 3.20453, 8.79773, 33.9392, 142.047, 633.503),
    # },
    "step06_vpbjacobi_dispscale": {
        "n_iterations": (12, 16, 17, 20, 25, 31),
        "total_time": (5.58084, 7.02474, 14.2053, 43.2277, 166.428, 716.287),
        "solve_time": (1.86589, 3.0295, 7.58369, 27.4059, 112.378, 498.27),
    },
}

LABELS = (
    "fieldsplit (selfp)",
    "vpbjacobi (defaults)",
    "vpbjacobi (tunefine)",
    "vpbjacobi (dispscale)",
)

pyplot.style.use("matplotlib_extras.color-lightbg")
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
