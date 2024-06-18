(sec-install)=
# Installation

{numref}`fig:install:choices` provides a guide to select the appropriate method for installing PyLith.
Installation of PyLith on a desktop or laptop machine is, in most cases, very easy.
Binary packages have been created for Linux and macOS platforms.
For Windows users, we recommend installing the Windows Subsystem for Linux and using the Linux binary (see instructions in Section {ref}`sec:install:windows`).
You can also run PyLith inside a Docker container, which provides a virtual Linux environment on any platform that Docker supports, including Linux, macOS, and Windows.
Installation of PyLith on other operating systems - or installation on a cluster - requires building the software from the source code, which can be difficult for inexperienced users.
We have created a small utility called PyLith Installer that makes installing PyLith and all of its dependencies from source much easier.

:::{figure-md} fig:install:choices
<img src="figs/install_choices.*" alt="Guide for selecting the appropriate installation choice based on a hardware and intended use. The installation options are discussed in more detail in the following sections." width = "100%" />

Guide for selecting the appropriate installation choice based on a hardware and intended use.
The installation options are discussed in more detail in the following sections.
:::

Help for installing and using PyLith is available from both a CIG forum and the GitHub issue tracking system <https://github.com/geodynamics/pylith/issues>.
See {ref}`sec-getting-help` for more information.

## Binary Package

The binaries are intended for users running on laptops or desktop computers (as opposed to clusters).
The binaries contain the compilers and header files, so users wishing to extend the code can still use the binary and do not need to build PyLith and its dependencies from source.
See {ref}`sec-developer-contributing` for more information on extending PyLith.

Binary executables are available for Linux (glibc 2.12 and later) and macOS (Intel 10.15 and later and arm64 11.0 and later) from the PyLith web page <https://geodynamics.org/resources/pylith/>.
Users running Windows 10 build 14316 and later can install a Linux bash environment and use the PyLith binary for Linux (see Section {ref}`sec:install:windows` for more information).

:::{tip}
On Linux systems you can check which version of glibc you have by running `ldd-version`

On macOS systems you can check the operating system version by clicking on the Apple icon and *About This Mac*.
:::

### Linux and macOS

1. Open a terminal window and change to the directory where you want to place the distribution.
    ```{code-block} bash
    $ cd $HOME
    $ mkdir pylith
    $ cd pylith
    ```
2. Download the Linux or macOS tarball from the PyLith web page
    <https://geodynamics.org/resources/pylith/supportingdocs/>, and save it to
    the desired location, e.g., `$HOME/pylith`.
3. Unpack the tarball.
    ```{code-block} bash
      # Linux 64-bit
      $ tar -xzf pylith-4.1.2-linux-x86_64.tar.gz

      # macOS
      $ tar -xzf pylith-4.1.2-macOS-10.15-x86_64.tar.gz
      ```
4. Set environment variables.
The provided `setup.sh` script only works if you are using bash shell.
If you are using a different shell, you will need to alter how the environment variables are set in `setup.sh`.
```{code-block} bash
$ source setup.sh
Ready to run PyLith.
```

:::{tip}
To bypass macOS quarantine restrictions, simply use command line program `curl` to download the tarball from within a terminal rather than using a web browser.

```{code-block} console
curl -L -O https://github.com/geodynamics/pylith/releases/download/v4.1.2/pylith-4.1.2-macOS-10.15-x86_64.tar.gz
```

Alternatively, if you do download the tarball using a web browser, after you unpack the tarball you can remove the macOS quarantine flags using the following commands (requires Administrator access):

```{code-block} console
# Show extended attributes
xattr ./pylith-4.1.2-macOS-10.15-x86_64

# Output should be
com.apple.quarantine

# Remove quarantine attributes
sudo xattr -r -d com.apple.quarantine ./pylith-4.1.2-macOS-10.15-x86_64
```
:::

:::{warning}
The binary distribution contains PyLith and all of its dependencies.
If you have any of this software already installed on your system, you need to be careful in setting up your environment so that preexisting software does not conflict with the PyLith binary.
By default the `setup.sh` script will prepend to the PATH and PYTHONPATH (for macOS and Linux) and LD_LIBRARY_PATH (for Linux) environment variables.
This will prevent most conflicts.
:::

(sec:install:windows)=
### Windows Subsystem for Linux Setup

PyLith is developed within the Unix/Linux framework, and we do not provide a native PyLith binary distribution for Windows.
The preferred approach for installing PyLith on a computer running Windows is to enable use of a Linux subsystem.
This permits use of the PyLith Linux x86_64 binary within the bash environment.

To enable the Linux subsystem on Windows 10 build 14316 and later (users running an earlier Windows build should use the PyLith Docker container):

1. Go to *Settings* &#8594; *Security* .
2. Under *For developers* select *Developer mode*.
This step should not be required for Windows build 16215 and later.
3. Go to *Control Panel* &#8594; *Programs* &#8594; *Turn Windows Features On or Off*.
4. Enable *Windows Subsystem for Linux* and click *OK* .
5. Restart the computer.
6. Go to *Start* &#8594; *bash*.
You will be prompted to download "Bash on Ubuntu on Windows" from the Windows Store.
Create a user account and password for the bash environment.
7. Install the PyLith Linux x86 binary within the bash environment following the instructions for installing the PyLith binary for Linux.
You will run PyLith within the bash environment just like you would for a Linux operating system.

#### WSL Setup for Gmsh and PyVista

The PyLith Linux binary package includes Gmsh and PyVista, which rely on several Linux libraries that are not included with PyLith.
In this section, we outline how to install these libraries based on user experiences; we do not have access to a Windows Subsystem for Linux for testing.

```{code-block} console
---
caption: Installation of extra libraries needed for Gmsh and PyVista. You only need to run these commands once per WSL installation.
---
# Update Linux
sudo apt-get update

# Install X libraries
sudo apt-get install libgl1-mesa-dri libglu1-mesa x11-apps gnome-terminal -y

# Make drivers visible to locations where Gmsh and PyVista expect them to be (/usr/lib/dri)
sudo ln -s /usr/lib/x86_64-linux-gnu/dri /usr/lib/dri

# Install libstdc++6
sudo apt-get install libstdc++6
```

:::{tip}
If you are able to run Gmsh or pylith_viz, but do not see any graphics, then you likely need to turn on software rendering by setting the `LIBGL_ALWAYS_SOFTWARE` environment variable.
You can set this environment variable using the `enable-software-rendering` argument to the PyLith `setup.sh` script or set it in your shell.

```{code-block} bash
# Turn on libGL software rendering.
export LIBGL_ALWAYS_SOFTWARE=1
```
:::

:::{seealso}
Refer to [Microsoft Tutorial: Run Linux GUI apps on the WSL](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gui-apps) for additional information about running graphic user interface applications within the Windows Subsystem for Linux.
:::

### Extending PyLith or Integrating Other Software Into PyLith

:::{note}
New in v3.0.0
:::

We strongly recommend using the [PyLith development environment Docker container](https://pylith-installer.readthedocs.io/en/latest/devenv/index.html) if you want to extend PyLith or integrate PyLith into other software.

## Source Installation

PyLith depends on a number of other packages (see {numref}`fig:pylith:dependencies`).
This complicates building the software from the source code.
In many cases some of the packages required by PyLith are available as binary packages.
On the one hand, using the binary packages for the dependencies removes the burden of configuring, building, and installing these dependencies, but that can come with its own host of complications if consistent compiler and configuration settings are not used across all of the packages on which PyLith depends.
This is usually not an issue with Linux distributions, such as Fedora, Ubuntu, and Debian that have good quality control; it can be an issue with macOS package managers, such as Fink, MacPorts, and Homebrew, where there is limited enforcement of consistency across packages.
Nevertheless, PyLith can be built on most systems provided the instructions are followed carefully.
PyLith is developed and tested on Linux and macOS.

A small utility, PyLith Installer, removes most of the obstacles in building PyLith and its dependencies from source.
For each package this utility downloads the source code, configures it, builds it, and installs it.
 This insures that the versions of the dependencies are consistent with PyLith and that the proper configure arguments are used.
 The minimum requirements for using the PyLith installer are a C compiler, `tar`, and `wget` or `curl`. Detailed instructions for how to install PyLith using the installer are included in the installer distribution, which is available from the PyLith web page <https://geodynamics.org/resources/pylith/supportingdocs/>.

## Verifying PyLith Installation

The easiest way to verify that PyLith has been installed correctly is to run one or more of the examples supplied with the binary and source code.
In the binary distribution, the examples are located in `src/pylith-4.1.2/examples` while in the source distribution, they are located in `pylith-4.1.2/examples`.
{ref}`sec-examples` discusses how to run and visualize the results for the examples.
To run the example discussed in Section {ref}`sec-examples-box-2d`:

```{code-block} bash
$ cd examples/box-2d
$ pylith step01_axialdisp.cfg
# A bunch of stuff will be written to stdout. The last few lines should be:
 >> .../lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

If you run PyLith in a directory without any input, you will get the error
message:

```{code-block} bash
$ pylith
 >> {default}::
 -- pyre.inventory(error)
 -- metadata.description <- ''
 -- Nonempty string required.
 >> {default}::
 -- pyre.inventory(error)
 -- metadata.arguments <- '[]'
 -- List of command line arguments required.
 >> {default}::
 -- pyre.inventory(error)
 -- metadata.pylith_version <- '[]'
 -- List of PyLith version constraints required.
 >> {default}::
 -- pyre.inventory(error)
 -- meshimporter.meshioascii.filename <- ''
 -- Filename for ASCII input mesh not specified.  To test PyLith, run an example as discussed in the manual.
 >> {default}::
 -- pyre.inventory(error)
 -- timedependent.problem_defaults.name <- ''
 -- Missing required property 'name' in default options for problem.
pylithapp: configuration error(s)
```

This indicates that at a very minimum metadata and the finite-element mesh file must be specified in order to run PyLith.

## Configuration on a Cluster

If you are installing PyLith on a cluster with a batch system, you can configure Pyre such that the `pylith` command automatically submits jobs to the batch queue.
Pyre contains support for the LSF, PBS, SGE, and Globus batch systems.
Properly configured, Pyre can handle job submissions automatically, insulating users from the details of the batch system and the site configuration.

See {ref}`sec-run-pylith-cluster` for more information on how to use PyLith on a cluster with a job scheduling system.
