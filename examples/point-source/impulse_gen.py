import numpy as np

def generate_ramp_stf(duration=0.5, dt=0.01, max_amplitude=1.0):
    time = np.arange(0, duration + dt, dt)
    amplitude = max_amplitude - (max_amplitude / duration) * time
    return time, amplitude

# Parameters
duration = 0.10
dt = 0.001
max_amp = 1.0
time, amplitude = generate_ramp_stf(duration, dt, max_amp)

# Write in the same format as impulse.timedb
with open("ramp.timedb", "w") as f:
    f.write("// -*- C++ -*- (tell Emacs to use C++ mode for syntax highlighting)\n")
    f.write("//\n")
    f.write("// Ramp source time function with linear increase to max amplitude.\n")
    f.write("#TIME HISTORY ascii\n")
    f.write("TimeHistory {\n")
    f.write(f"  num-points = {len(time)+1} // number of points in time history\n")
    f.write("  time-units = second // units for time\n")
    f.write("}\n")
    for t, a in zip(time, amplitude):
        f.write(f"{t:10.4f}   {a:10.4f}\n")
    t_final = 999.000
    a_final = 0.000
    f.write(f"{t_final:10.4f}   {a_final:10.4f}\n")

