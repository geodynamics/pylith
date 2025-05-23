(sec-developer-workflow)=
# Overview

We use the Git version control system <https://git-scm.com> with a central GitHub repository <https://github.com/geodynamics/pylith> for development.
We will refer to this central repository as the `geodynamics/pylith` repository.
Only the PyLith maintainers have write access to the `geodynamics/pylith` repository; everyone else is limited to read access.
Contributions from the community are incorporated into the `geodynamics/pylith` repository via pull requests.

Currently, the PyLith maintainers use the
[Forking Workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow).
Each developer forks the geodynamics repository and a new feature is added on its own branch in the developer fork.
Once it is implemented, tested, and documented, it is merged to the `main` branch of the `geodynamics/pylith` repository.

The forking workflow  is illustrated in Figures [git-repositories](fig-developer-git-repositories) and [git-branch](fig-developer-git-branch).
You fork the `geodynamics/pylith` repository so that you have a copy of the repository for your changes.
Each new, independent feature you develop should be done in a separate branch.
Large feature contributions should be broken up into multiple phases, where each phase provides a working version that passes all of its tests.
Each phase corresponds to a separate branch that is merged back into the `geodynamics/pylith` repository (usually the `main` branch) via a [pull request](https://help.github.com/articles/about-pull-requests/) when it is completed.

:::{figure-md} fig-developer-git-repositories
<img src="figs/gitworkflow_repositories.*" alt="" width="80%" />

Overview of repositories for the Git forking workflow used in PyLith development.
The main repository is `geodynamics/pylith` at GitHub.
Developers create a fork of that repository under their own GitHub account (e.g., `saradeveloper`), which is cloned to their local computer.
Changes and additions to the code are committed to the repository on the local computer, which are pushed to the developer's GitHub account.
Once a development task is completed, a developer is encouraged to contribute the changes and additions to the main repository via pull requests.
:::

:::{figure-md} fig-developer-git-branch
<img src="figs/gitworkflow_branch.*" alt="" width="100%" />

Overview of repositories and branches for the Git forking workflow used in PyLith development.
You keep the `main` branch on your local machine in sync with the `main` branch in the community repository; we never use the `main` branch in your GitHub repository.
From the `main` branch on your local machine, you create a feature branch, e.g., `feature-powerlaw`, to complete a single task such as adding a feature, fixing a bug, or making an improvement to the documentation.
You should break up the changes into several steps, which are saved locally as commits.
The commits are pushed to your repository on GitHub for backup or syncing across multiple computers.
Once the task is completed, you submit a pull request to have the changes merged to the `main` branch in the community repository.
Once the pull request is merged, you update your local `main` branch from the community repository and then push the changes to your repository on GitHub.
:::
