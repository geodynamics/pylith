# Git Quick Reference

## Break commit into multiple commits

```{code-block} bash
git rebase -i main
# Set commit to 'edit', exit, and save.
git reset HEAD~
# Make changes.
git rebase --continue
```

## Update from upstream

1. Set remote (one time)

```{code-block} bash
# Show current remotes
git remote -v

# Add remote upstream
git remote add upstream UPSTREAM_REPO

# Verify addition of upstream remote
git remote -v
```

2. Update upstream repo

```{code-block} bash
git fetch upstream -p
```

## Rebase branch off `main`

```{code-block} bash
# Make sure main is current
git pull main

# Check out branch you want to rebase
git checkout BRANCH

# Rebase branch with respect to `main`
git rebase -i main
```

## Branches

* **Checkout branch and set tracking**: `git checkout -b $BRANCH --track $REMOTE/$BRANCH`
* **Compare current heads**: `git diff branch1..branch2`
* **Compare with respect to common ancestor**: `git diff main...feature`
* **Compare commits**: `git log branch1..branch2`
* **Compare specific file**: `git diff main..feature -- FILE`
* **Compare specific file wih working tree**: `git diff main-- FILE`

### Delete branches

* **Local**: `git branch -d BRANCH`
* **Remote**: `git push origin --delete BRANCH`

### Rename local branch

`git branch -m OLD NEW`

### Show merged/unmerged branches

```{code-block} bash
git branch --merged $BRANCH
git branch --no-merged $BRANCH
```
