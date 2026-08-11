# Releasing pyintval

This project ships prebuilt wheels to PyPI via GitHub Actions using PyPI
**Trusted Publishing** (OIDC), so no API tokens are stored anywhere. This guide
covers the one-time setup and the per-release process.

## One-time setup

### 1. Create the GitHub repository and push

```sh
# from the project root
git remote add origin https://github.com/marciogameiro/pyintval.git
git branch -M main
git push -u origin main
```

CI (`.github/workflows/ci.yml`) runs on every push and pull request. The
release workflow (`.github/workflows/release.yml`) also builds wheels on pull
requests that touch packaging files, so you can confirm the ARM64 and Windows
(clang-cl) builds are green **before** tagging a release.

### 2. Configure the PyPI trusted publisher

On <https://pypi.org> (log in, then go to your account → *Publishing*, or the
project's *Settings → Publishing* once it exists), add a **pending publisher**
with:

| Field | Value |
|---|---|
| PyPI project name | `pyintval` |
| Owner | `marciogameiro` |
| Repository name | `pyintval` |
| Workflow name | `release.yml` |
| Environment name | `pypi` |

Repeat on <https://test.pypi.org> with environment name `testpypi` if you want
TestPyPI dry-runs.

### 3. Create the GitHub environments

In the repo: *Settings → Environments → New environment*, create `pypi` (and
optionally `testpypi`). These match the `environment:` names in the release
workflow and are where you can add protection rules (e.g. require a reviewer
before publishing).

### 4. Enable Read the Docs

Import the repo on <https://readthedocs.org>; it reads `.readthedocs.yaml`, which
**builds the package from source** with GCC (so `autodoc` documents the exact
code of the branch or tag being built — no dependency on the published wheel and
no release-time race) and publishes the Sphinx docs.

Then configure the project so releases publish on their own:

- **Versions** tab → activate the **`stable`** version. It auto-tracks the
  newest SemVer tag and re-points itself whenever a newer tag is built.
- **⚙ Admin → Settings → Default version** → **`stable`**, so the docs home page
  shows the newest *release* (`latest` keeps tracking `main` for dev docs).
- **⚙ Admin → Automation Rules → Add rule:** match *SemVer versions*, version
  type **Tag**, action **Activate version** — so every new `vX.Y.Z` tag is
  activated and built automatically.

The docs badge/URL in `README.md` already point at the RTD site.

## Cutting a release

1. Bump the version in `pyproject.toml`, and move the `## [Unreleased]` section
   of `CHANGELOG.md` to the new version (add the date and a release link). If the
   `Development Status` classifier is still `3 - Alpha`, consider bumping it to
   `4 - Beta`.
2. Commit: `git commit -am "Release X.Y.Z"`.
3. (Optional) Dry-run to TestPyPI: from the Actions tab, run the **Release**
   workflow manually (`workflow_dispatch`) with `publish_target = testpypi`,
   then `pip install -i https://test.pypi.org/simple/ pyintval` in a fresh venv.
4. Tag and push:
   ```sh
   git tag v0.1.0
   git push origin v0.1.0
   ```
   The tag triggers the full matrix (Linux x86-64/aarch64, macOS arm64/x86-64,
   Windows AMD64 via clang-cl, plus the sdist), tests every wheel against the
   suite, and publishes to PyPI.
5. Create a GitHub Release from the tag with notes (optional but recommended).

Documentation updates itself: the tag push makes Read the Docs build that version
from source and re-point `stable` (the default version) to it, so the docs home
page shows the new release with no manual step. Merging to `main` separately
rebuilds `latest` (the dev docs). See *One-time setup → Enable Read the Docs*.

## Notes

- Bump the `pypa/cibuildwheel@vX.Y.Z` pin in `release.yml` periodically to pick
  up support for new Python versions.
- The Windows wheels **must** be built with clang-cl; the workflow selects it
  via `CMAKE_GENERATOR_TOOLSET=ClangCL`. A plain-MSVC build fails fast with an
  explanatory CMake error.
- CPython 3.14 wheels are built once the pinned cibuildwheel version supports
  that tag.
