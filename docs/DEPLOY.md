# Deploying documentation to https://ciera-transients.github.io/POTPyRI/

The docs are built with Sphinx and deployed to GitHub Pages via the `gh-pages` branch. The URL is a **project site**: `https://<org>.github.io/<repo>/` → `https://ciera-transients.github.io/POTPyRI/`.

## One-time setup (repo maintainers)

1. **Enable GitHub Pages**
   - On GitHub: **Settings** → **Pages** (left sidebar).
   - Under **Build and deployment**:
     - **Source**: Deploy from a branch.
     - **Branch**: choose `gh-pages`, folder **/ (root)**.
     - Save.

2. **First deployment**
   - The workflow [`.github/workflows/documentation.yml`](../.github/workflows/documentation.yml) runs on every **push to `main`**.
   - It builds the docs (`pip install .[docs]`, `cd docs && make html`) and pushes `docs/build/html/` to the `gh-pages` branch.
   - After the first successful run, the site is available at https://ciera-transients.github.io/POTPyRI/ (can take a few minutes).

3. **Pull requests**
   - The same workflow runs on pull requests targeting `main`, but the **Deploy** step is skipped (deploy only runs on push to `main`). So PRs validate that the docs build without updating the live site.

## Local build

To build and view the docs locally:

```bash
pip install .[docs]
cd docs && make html
open build/html/index.html   # or open build/html/index.html on Linux
```

## Configuration

- **Sphinx config**: `docs/source/conf.py` (e.g. `html_baseurl`, theme, extensions).
- **Base URL**: `html_baseurl = 'https://CIERA-Transients.github.io/POTPyRI/'` is already set so links and assets resolve correctly on GitHub Pages.

If the site returns 404 after enabling Pages, wait a few minutes and ensure the `gh-pages` branch exists and has content (run the workflow once by pushing to `main`).
