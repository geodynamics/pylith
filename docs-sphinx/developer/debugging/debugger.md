# Using the debugger

Please see the `gdb` and `lldb` documentation for detailed instructions.
Here we illustrate some common basic commands.

```{code-block} console
---
caption: Debugging with gdb
---
# Set breakpoint at line 150 of Material.cc
(gdb) b Material.cc:150

# Set breakpoint at exception throw
(gdb) catch throw

# Show arguments for the current frame
(gdb) info args

# Show local variables for the current frame
(gdb) info locals

# Show the contents of a local variable: p VARIABLE
(gdb) p numFields

# Show the contents of local array: p POINTER[0]@SIZE
# Print array of 4 values pointed to by variable values
(gdb) p values[0]@4

# Print stack trace
(gdb) backtrace
```

```{code-block} console
---
caption: Debugging with lldb
---
(lldb) b Material.cc:150

# Set breakpoint at exception throw
(lldb) break set -E C++

# Show local variables
(lldb) frame variable

# Show the contents of a local variable: frame variable VARIABLE
(lldb) frame variable numFields
# Alternatively
(lldb) p numFields

# Show the contents of an array of values: parray SIZE POINTER
# Show the contents of array of 10 values pointed to by the variable values.
(lldb) parray 10 values
```

## C++ and MMS tests

The executables in the build directory are shell script wrappers created by `libtool`.
The underlying binary executables are in the `.libs` directory.
When using the debugger, pass the binary executable to the debugger.
For example, `gdb .libs/test_problems`.

## Full-scale tests

To start the `gdb` debugger when running the PyLith application, simply add the command line argument `--petsc.start_in_debugger`.
To use an alternative debugger, such as `lldb`, append the name of the debugger executable, for example `--petsc.start_in_debugger=lldb`.
By default, PETSc will try to start the debugger in an xterm.
To use an alternative terminal program, use the command line argument `--petsc.debug_terminal=TERMINAL`.
For example for the GNOME terminal, use `--petsc.debug_terminal="gnome-terminal -x"`.
