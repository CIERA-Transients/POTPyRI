#!/bin/bash
# Trim POTPyRI git history by removing large files/dirs from all commits.
# Run from repo root. Back up the repo first (e.g. clone to another folder).
# After running: git push --force-with-lease (all branches must be force-pushed).

set -e
cd "$(git rev-parse --show-toplevel)"

# Paths to remove from entire history:
# - Outside potpyri/: old calibrations/, specpipe/, staticmasks/ (not in current tree)
# - Obsolete under potpyri/: potpyri/data/calibrations/ (replaced by potpyri/data/cal/)
PATHS_TO_REMOVE=(
  "calibrations"
  "specpipe"
  "staticmasks"
  "potpyri/data/calibrations"
)

# Install git-filter-repo if missing
if ! command -v git-filter-repo &>/dev/null; then
  echo "Installing git-filter-repo..."
  pip install git-filter-repo
fi

# Build --path args (all paths then single --invert-paths)
PATH_ARGS=()
for p in "${PATHS_TO_REMOVE[@]}"; do
  PATH_ARGS+=(--path "$p")
done
PATH_ARGS+=(--invert-paths)

echo "Removing paths from history: ${PATHS_TO_REMOVE[*]}"
# Non-interactive: set TRIM_HISTORY_CONFIRM=y to skip prompt
if [[ "${TRIM_HISTORY_CONFIRM}" != "y" ]]; then
  echo "This rewrites history. Ensure you have a backup."
  read -p "Continue? [y/N] " -n 1 -r
  echo
  [[ $REPLY =~ ^[yY]$ ]] || exit 1
fi

git filter-repo "${PATH_ARGS[@]}" --force

echo "Running git gc --aggressive..."
git reflog expire --expire=now --all
git gc --prune=now --aggressive

echo "Done. Next: git remote add origin <url>  # if remote was removed"
echo "Then: git push --force-with-lease origin main dev-cdk  # and other branches"
