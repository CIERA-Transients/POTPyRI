# Git history trim (large files removed)

History was rewritten to remove the following paths from **all commits**:

- **`calibrations/`** (root) – old BINOSPEC/MMIRS calibration FITS
- **`specpipe/`** (root) – legacy spectral reduction tree and LRIS cals
- **`staticmasks/`** (root) – old large static mask FITS (replaced by `potpyri/data/staticmasks/`)
- **`potpyri/data/calibrations/`** – old layout (replaced by `potpyri/data/cal/`)

**Result:** `.git` size went from ~785 MB to ~569 MB.

## After pulling this repo (or if you ran the trim locally)

1. **Do not run `git fetch origin` yet** – it would re-download the old history and undo the size reduction. Remote-tracking refs were removed so that pruned objects could be dropped.

2. **Force-push rewritten branches** (history changed; everyone must use the new history):
   ```bash
   git push --force-with-lease origin main
   git push --force-with-lease origin dev-cdk
   ```
   Push any other branches you care about the same way.

3. **Then** restore remote tracking refs and get any updates:
   ```bash
   git fetch origin
   ```

4. **Tags** were rewritten as well. To update them on the remote:
   ```bash
   git push origin --force --tags
   ```

6. **Collaborators** should re-clone or run:
   ```bash
   git fetch origin
   git reset --hard origin/<branch>
   ```
   for each branch they use (old history is incompatible).

## Running the trim again (e.g. on another clone)

See `scripts/trim_history.sh`. You can run it with:

```bash
TRIM_HISTORY_CONFIRM=y ./scripts/trim_history.sh
```

Or use `git filter-branch` as was done (paths are listed in the script).
