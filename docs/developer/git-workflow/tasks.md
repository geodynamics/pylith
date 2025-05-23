# Git tasks

## Setting up your fork

There are two steps to setting up a copy of the PyLith repository that you can use for development under the Forking Workflow:

1. Create a GitHub account.
2. Fork the `geodynamics/pylith` repository.

    This will create a copy of the PyLith repository in your GitHub account.
    You can make changes to this copy and, when you are ready to contribute changes back to the `geodynamics/pylith` repository, you will create a pull request.

:::{important}
Create a personal access token for your GitHub account.
See [GitHub Docs: Creating a personal access token](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token) for more information.
This provides additional security and eliminates the need for you to enter your username and password whenever you push to your repository.

We strongly recommend that you setup your GitHub account to use two-factor authentication.
See [GitHub Docs: Securing your account with two-factor authentication](https://help.github.com/articles/securing-your-account-with-two-factor-authentication-2fa/).

You should also sign your commits using GPG or S/MIME.
See [GitHub Docs: Managing commit signature verification](https://docs.github.com/en/github/authenticating-to-github/managing-commit-signature-verification) for more information.
:::

### Set the upstream repository

For each clone of your fork (each computer with a local copy of your fork), you need to create a link to the "upstream" `geodynamics/pylith` repository.
This allows you to keep your repository in sync with the community repository.
If you are using the PyLith Docker development container and followed all of the installation instructions, then you already completed this step; you can list the current remotes for your fork to verify.

```{code-block} bash
---
caption: Setting upstream repository
---
# List the current remotes for your fork.
git remote -v
# Output
origin https://github.com/YOUR_GITHUB_USERNAME/pylith.git (fetch)
origin https://github.com/YOUR_GITHUB_USERNAME/pylith.git (push)

# Set the link to the remote upstream repository
git remote add upstream https://github.com/geodynamics/pylith.git

# Verify the upstream repository has been added.
git remote -v
# Output:
origin https://github.com/YOUR_GITHUB_USERNAME/pylith.git (fetch)
origin https://github.com/YOUR_GITHUB_USERNAME/pylith.git (push)
upstream https://github.com/geodynamics/pylith.git (fetch)
upstream https://github.com/geodynamics/pylith.git (push)
```

(sec-developer-feature-branch)=
## Creating a new feature branch

```{code-block} bash
---
caption: Creating a feature branch
---
# Start from the current development branch (usually "main")
git checkout main

# Make sure it is up to date.
git pull

# Create a new branch from 'main'.
git checkout -b BRANCH-NAME

# Examples
git checkout -b feature-powerlaw-rheology
git checkout -b fix-fault-output
```

:::{tip}
If you are implementing a feature that requires a significant amount of code development, we strongly recommend that you break the implementation into pieces that can each be tested, documented, and integrated into the PyLith `main` branch.
Another approach (equally valid) is to create a series of feature branches implementing the different phases that all get merged into the main feature branch; the main feature branch would then be merged into the `main` branch via a pull request.
:::

:::{tip}
Feature branches are a great way to experiment with an implementation.
You can create a feature branch and if you decide the implementation is headed in the wrong direction, you can simply create a new feature branch from your original starting point and delete the bad feature branch.
:::

## Staging, Committing, and Pushing Changes

The Git `add` and `commit` commands are used to stage and commit changes to a branch.
Staging refers to assembling the list of files to include in a commit.
A commit adds code changes to the current branch along with a message describing the changes.
A commit changes only the current branch on your local machine.
In order to update your GitHub repository you need to `push` your changes.
See the Git documentation for details about these commands.
There are Git interfaces built into a number of editors and integrated development environments (IDEs), as well as standalone graphical user interfaces to Git.

:::{tip}
Commit messages should explain changes.
Keep the first line under 80 characters, if possible, and include an empty line between the first line and any additional lines that provide more detailed explanations.

The commit messages provide important documentation on why code is changed.
Your future self and your fellow developers will appreciate good explanations of all changes.
:::

:::{warning}
If you have multiple branches, make sure you are on the correct branch before making commits.
:::

:::{tip}
You should only make commits on your local feature branches.
These are also the only branches that you should push to your GitHub repository.

We strongly recommend never pushing branches from the upstream repository to your GitHub repository.
:::

## Keeping a feature branch in sync with branch `main`

The PyLith `main` branch may change while you are working on a feature branch.
In some cases it might have a new feature that you do not need and does not affect your work, in which case it does not matter whether you incorporate those changes into your feature branch.
However, in other cases, there might be an important bugfix or feature that you want to use in your feature branch while you are working on it.
The recommended way to incorporate these changes is to rebase your feature branch relative to the PyLith `main` branch.

Rebasing essentially replays your commits on top of the commits in the other branch.
With interactive rebasing you can also rewrite the commit history of a feature branch by reordering, dropping, and/pr combining (squashing) commits.
See [Git: Rewriting History](https://git-scm.com/book/en/v2/Git-Tools-Rewriting-History) for more information.

:::{important}
Make sure all of your local changes have been committed or [stashed](https://git-scm.com/docs/git-stash).
:::

:::{danger}
Rebasing and force pushing can cause irreversible damage to your repository.
We recommend practicing rebasing using a toy repository before attempting to rebase with real code.
:::

The primary steps to rebasing with respect to the PyLith `main` branch are:

1. Commit or stash any changes so that you do not have any modified files.
2. Start the interactive rebasing by running `git rebase -i main`.
3. Edit the commit history as desired.
4. Fix any conflicts during the rebase.
5. Verify the rebase is correct by running tests, etc.
6. Force push the branch using `git push --force`.

:::{danger}
After you have rebased but before doing a forced push, running `git status` will report something like:

```bash
On branch saradeveloper/my-great-new-feature
Your branch and 'origin/saradeveloper/my-great-new-feature' have diverged,
and have 56 and 16 different commits each, respectively.
  (use "git pull" to merge the remote branch into yours)
```

Do not run `git pull`!!!!!
This is one of the few times that you should *not* do what git suggests.
Doing so will try to merge your original history into your new one.
Instead, you should verify that the new commit history is correct and then do a force push by running `git push --force`.
:::

:::{tip}
If anytime *during* the rebasing process, you make a mistake or decide you want to abort the rebasing process, simple run `git rebase --abort`.
This will return the repository to the state immediately before the rebasing process.
:::

(sec-developer-pull-request)=
## Creating a pull request

Once you have completed implementing, testing, and documenting a new feature on a branch in your fork and wish to contribute it back to the `geodynamics/pylith` repository, you open a pull request.
See [GitHub Docs: About pull requests](https://help.github.com/articles/about-pull-requests/) for more information.
When you create a pull request, GitHub will run a series of tests using a temporary copy with the merged changes using GitHub actions.
The GitHub webpage for the pull request will have a `Checks` tab showing the progress of the tests and any errors.
One or more of the PyLith maintainers will review your pull request and may request changes.
You can continue to make commits to your branch; the pull request will be automatically updated whenever you push the changes to your GitHub repository.

:::{tip}
Be sure the checkbox `Allow edits and access to secrets by maintainers` is checked so maintainers can make corrections.
:::

:::{tip}
To become familiar with making pull requests, we recommend starting with a small, simple change.
This may be as little as fixing a typo in the documentation or a comment.
Create a feature branch for the change, push it to your repository, and then make a pull request.
:::

## Adding remotes for accessing other PyLith forks

When collaborating with other people working on PyLith, it is helpful to be able to checkout branches from their forks.
You can add their fork as an additional "remote" repository.

```{code-block} bash
---
caption: Adding an additional remote repository to track branches in other forked repositories.
---
# Add remote
git remote add FORK-NAME https://github.com/GITHUB_USERNAME/pylith.git
# Example:
git remote add saradeveloper https://github.com/saradeveloper/pylith.git

# Show remotes
git remote -v

# Fetch the information for the remote
git fetch FORK-NAME
# Example:
git fetch saradeveloper

# Checkout remote branch
git checkout -b saradeveloper/feature-powerlaw-rheology

# Push to remote branch (requires write access)
git push FORK-NAME BRANCH
# Example:
git push saradeveloper feature-powerlaw-rheology
```
