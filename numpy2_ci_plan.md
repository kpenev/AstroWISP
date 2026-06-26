# Adding an OS × Python × numpy-1/2 test matrix to AstroWISP CI

## Goal

Today `.github/workflows/build_wheels.yml` only *builds* wheels (one job,
`workflow_dispatch`-only, currently restricted to `ubuntu-latest`). It never
imports the result or runs a single test. Before we can lift the `numpy < 2`
pin and let AutoWISP move to Python 3.13, we need CI that proves the compiled
extension works across:

* **OS**: `ubuntu-latest`, `windows-2022`, `macos-13` (Intel), `macos-14` (Apple silicon)
* **Python**: 3.11, 3.12, 3.13, 3.14
* **numpy**: a 1.x line (`numpy<2`) and a 2.x line (`numpy>=2`)

## Facts that shape the design

1. **The wheel is numpy-agnostic.** AstroWISP talks to its C++ library through
   `ctypes` (`numpy.ctypeslib.ndpointer`), *not* the numpy C-API. numpy is not a
   build dependency and is not baked into the `.so`/`.pyd`. So one wheel runs
   under numpy 1 *or* 2 — the numpy axis only changes the **test** environment,
   never the build. (This is why the matrix below varies numpy in
   `CIBW_TEST_REQUIRES`, not in the build.)

2. **Tests and their data ship in the wheel.** `PythonPackage/astrowisp/tests/`
   and `tests/test_data/` are installed via `py.install_sources` / `subdir(...)`
   in meson. So cibuildwheel can run the suite against the *installed* package
   with `pytest --pyargs astrowisp.tests` — no source checkout needed in the
   test venv, and no test-data download. (`test_grcollect.py` is currently kept
   out of `tests/meson.build` so its hard-coded `B:/...` path never gets
   collected; the "drop the grcollect executable and its test" prerequisite
   below removes it outright.)

3. **Boost is the only heavy native dep**, and it is already handled per-OS in
   the existing job (MarkusJx/install-boost for mac/windows,
   `CIBW_BEFORE_BUILD_LINUX` yum for the manylinux container). cibuildwheel's
   built-in test support reuses all of that, so we should extend the existing
   job rather than write a from-source test job that re-solves Boost.

4. **numpy/Python compatibility is not a full rectangle** — numpy 1.26 (the last
   1.x) has no cp313 wheels. So the (numpy 1 × Python 3.13) cell is impossible
   and must be excluded.

   | Python | numpy 1.26 | numpy ≥ 2 |
   |--------|:----------:|:---------:|
   | 3.11   | ✅          | ✅         |
   | 3.12   | ✅          | ✅         |
   | 3.13   | ❌ (no wheel)| ✅         |
   | 3.14   | ❌ (no wheel)| ✅         |

   The numpy-1 line tops out at 3.12; only the numpy-2 line reaches the newest
   interpreters (3.13, 3.14). **Caveat:** building cp313/cp314 requires a
   recent cibuildwheel — the workflow currently pins `cibuildwheel@v2.18.0`
   (mid-2024), which predates cp314. Bump it to cibuildwheel 3.x and confirm
   the C++ build + Boost-install action run on a 3.14 runner. AstroWISP loads
   its library via `ctypes`, so the compiled artifact isn't tied to the CPython
   C-API; the risk is in the build toolchain, not the runtime binding.

## Prerequisite (must land first)

`pyproject.toml` currently pins the runtime dependency to `numpy < 2`
(line 43). While that pin is in place, the numpy-2 cells cannot resolve
(`CIBW_TEST_REQUIRES: numpy>=2` conflicts with the package's own `numpy<2`).
Relax it to an unbounded `numpy` (keep a sane floor, e.g. `numpy>=1.21`) **in
the same PR** that adds this matrix. Bump `requires-python` from `>=3.7` to
`>=3.11` at the same time if AstroWISP is to match AutoWISP's floor.

## Prerequisite: drop the grcollect executable and its test

The `grcollect` executable is **not a viable tool**: its build fails on
Apple-silicon Macs, and it misbehaves on Windows. We therefore do not want to
build or test it on any platform — remove it entirely rather than ship a target
that is broken on half the matrix:

1. **Stop building/installing it.** Delete the second `executable(...)` block in
   `fitsh-0.9.4-minimum/src/meson.build` (the `'grcollect'` target, currently
   lines ~69-78) and the now-unused `grcollect_src` list (lines ~51-55). Only
   the `fistar` executable should remain installed into `astrowisp/`. The
   `grcollect.c` / `cache.c` source files can stay on disk; they just stop being
   compiled.

2. **Delete the test.** Remove
   `PythonPackage/astrowisp/tests/test_grcollect.py`. It is already absent from
   `PythonPackage/astrowisp/tests/meson.build`'s `py_src` list, so it is never
   shipped or collected — but the file (with its hard-coded
   `B:/fitsh-0.9.4-minimum/builddir/src/grcollect.exe` path) should go so nobody
   resurrects it. No change to `tests/meson.build` is needed since it was never
   listed there.

3. **Why now.** Lifting the numpy pin brings Apple-silicon (`macos-14`) and
   Windows into a green-by-default matrix, and `grcollect` cannot pass on either
   (broken build on arm64 Macs, wrong behavior on Windows). Dropping the target
   and its test removes a tool we've decided not to support, so the matrix
   reflects only what AstroWISP actually ships (`fistar` + the Python package).
   This also supersedes the parenthetical note in fact #2 above: rather than
   relying on `test_grcollect.py` being "intentionally not in `tests/meson.build`,"
   we delete it outright.

## Matrix mechanics (fused job — superseded, see "Decision" below)

cibuildwheel already iterates Python versions internally (one wheel per
`cp3XX`), and it can build a venv per wheel, install the wheel plus
`CIBW_TEST_REQUIRES`, and run `CIBW_TEST_COMMAND`. So we only add **two** matrix
axes' worth of behaviour:

* The **OS** axis becomes the full four-runner list (uncomment what's already
  there).
* A new **numpy** axis (`"1"`, `"2"`) drives both *which* Pythons cibuildwheel
  builds (`CIBW_BUILD`) and *which* numpy the test venv installs
  (`CIBW_TEST_REQUIRES`).

Net grid: `4 OS × {numpy1: cp311,cp312} ∪ {numpy2: cp311,cp312,cp313,cp314}` =
**24 build+test cells**, matching the table above with the impossible cells
(numpy 1 × 3.13/3.14) already excluded.

Add `fail-fast: false` so one red cell doesn't cancel the rest. Keep the
`workflow_dispatch`-only trigger (a deliberate choice in `build_wheels.yml` —
this grid is heavy and meant to be run on demand, not on every push).

> **Superseded:** the fused single-job approach described in this section
> rebuilds each wheel once per numpy axis (8 redundant builds). It is kept here
> only to explain the matrix mechanics. The **chosen** design is the two-job
> split in "Decision: split build (16 wheels) from test (24 environments)"
> below, which builds each wheel exactly once.

### Concrete workflow

Either fold this into `build_wheels.yml` or add a sibling
`.github/workflows/test_matrix.yml`. Key additions over the current file are
marked `# NEW` / `# CHANGED`:

```yaml
name: Test matrix

on: [workflow_dispatch]               # kept manual-only, as in build_wheels.yml

jobs:
  build_and_test:
    name: ${{ matrix.os }} · numpy ${{ matrix.numpy }}
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false                # NEW: don't cancel siblings on one failure
      matrix:
        os: [ubuntu-latest, windows-2022, macos-13, macos-14]  # CHANGED: full grid
        numpy: ["1", "2"]             # NEW: numpy axis

    steps:
      - uses: actions/checkout@v4

      - name: Detect Boost Platform
        if: ${{ matrix.os == 'macos-13' }}
        run: echo "BOOST_PLATFORM=10.15" >> "$GITHUB_ENV"

      - name: Detect Boost Platform
        if: ${{ matrix.os == 'macos-14' }}
        run: echo "BOOST_PLATFORM=14" >> "$GITHUB_ENV"

      - name: Install Boost Mac/Windows
        if: ${{ startsWith(matrix.os, 'mac') || startsWith(matrix.os, 'windows') }}
        id: install-boost
        uses: MarkusJx/install-boost@v2
        with:
            boost_version: ${{ matrix.os == 'macos-13' && '1.79.0' || '1.81.0' }}
            platform_version: ${{ startsWith(matrix.os, 'windows') && '2022' || env.BOOST_PLATFORM }}
            toolset: ${{ startsWith(matrix.os, 'mac') && 'clang' || 'mingw' }}
            arch: ${{ matrix.os == 'macos-14' && 'aarch64' || 'x86' }}
            boost_install_dir: ${{github.workspace}}/external

      - name: Build & test wheels
        uses: pypa/cibuildwheel@v3.0.0   # CHANGED: v2.18.0 predates cp313/cp314
        env:
            MACOSX_DEPLOYMENT_TARGET: ${{ matrix.os == 'macos-14' && '14.0' || '10.15' }}
            BOOST_ROOT: ${{ steps.install-boost.outputs.BOOST_ROOT }}

            # NEW: pick Pythons per numpy axis (numpy 1.26 stops at cp312)
            CIBW_BUILD: ${{ matrix.numpy == '1' && 'cp311-* cp312-*' || 'cp311-* cp312-* cp313-* cp314-*' }}

            # NEW: install the requested numpy line + pytest into the test venv,
            #      then run the installed suite (import-based, so cwd-independent)
            CIBW_TEST_REQUIRES: ${{ matrix.numpy == '1' && 'pytest "numpy<2"' || 'pytest "numpy>=2"' }}
            CIBW_TEST_COMMAND: python -m pytest --pyargs astrowisp.tests

            # --- unchanged build plumbing ---
            CIBW_BEFORE_BUILD_LINUX: >
              sed -i s/mirror.centos.org/vault.centos.org/g /etc/yum.repos.d/CentOS-*.repo &&
              sed -i s/^#.*baseurl=http/baseurl=http/g /etc/yum.repos.d/CentOS-*.repo &&
              sed -i s/^mirrorlist=http/#mirrorlist=http/g /etc/yum.repos.d/CentOS-*.repo &&
              yum -y -v install boost boost-thread boost-devel
            CIBW_BEFORE_BUILD_WINDOWS: pip install wheel delvewheel
            CIBW_REPAIR_WHEEL_COMMAND_MACOS: >
              DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$BOOST_ROOT/lib delocate-listdeps {wheel} &&
              DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$BOOST_ROOT/lib delocate-wheel --require-archs {delocate_archs} -w {dest_dir} {wheel}

      - name: Upload wheels
        if: ${{ matrix.numpy == '2' }}   # NEW: only keep one (numpy-agnostic) wheel set
        uses: actions/upload-artifact@v4
        with:
          name: cibw-wheels-${{ matrix.os }}
          path: ./wheelhouse/*.whl
```

### Why each piece

* **`CIBW_BUILD` keyed on `matrix.numpy`** is what excludes the impossible
  (numpy 1 × 3.13) cell without an `exclude:` block — under the numpy-1 axis we
  simply never build cp313.
* **`CIBW_TEST_REQUIRES` keyed on `matrix.numpy`** forces pip to install the
  chosen numpy line into cibuildwheel's per-wheel test venv. Because the wheel's
  own dep is just `numpy` (after the pin is relaxed), the constraint here decides
  1 vs 2.
* **`CIBW_TEST_COMMAND` uses `--pyargs astrowisp.tests`** so it tests the
  *installed* package (fact #2). Adjust if AstroWISP grows a `python -m
  astrowisp.tests` entry point — then `python -m astrowisp.tests` is equivalent
  and clearer.
* **`if: matrix.numpy == '2'` on upload** avoids uploading two identical wheel
  sets per OS (the wheel doesn't depend on numpy).

## Decision: split build (16 wheels) from test (24 environments)

**This is the chosen design** — the earlier fused "build once per numpy axis"
job is rejected. The wheel is numpy-agnostic, so building it once per numpy line
duplicates work: under the fused job cp311/cp312 get built twice (once for the
numpy-1 axis, once for numpy-2), 8 redundant builds across 4 OS. Instead build
each wheel **once** and fan it out to every applicable test environment. Two
jobs:

1. **`build` — 16 wheels.** One numpy-free build per (OS, Python), keeping
   separate arm64 and x86_64 mac wheels:

   | Build runner | Pythons | Wheels |
   |--------------|---------|:------:|
   | `ubuntu-latest` (manylinux x86_64) | cp311–cp314 | 4 |
   | `windows-2022` (amd64)             | cp311–cp314 | 4 |
   | `macos-13` (Intel x86_64)          | cp311–cp314 | 4 |
   | `macos-14` (Apple-silicon arm64)   | cp311–cp314 | 4 |
   | **Total**                          |             | **16** |

   This is the current cibuildwheel job run once per OS (no numpy axis),
   `upload-artifact` all 16. Reuses the existing per-OS Boost handling
   unchanged.

2. **`test` — 24 environments.** A plain `actions/setup-python` matrix over the
   four runner OSes × the six valid (Python, numpy) cells, `download-artifact`,
   `pip install <matching wheel> pytest asteval pandas "numpy<2|>=2"`, then
   `python -m pytest --pyargs astrowisp.tests`. Per runner OS the six cells are:

   | Python | numpy 1 | numpy 2 |
   |--------|:-------:|:-------:|
   | 3.11   | ✅       | ✅       |
   | 3.12   | ✅       | ✅       |
   | 3.13   | ❌       | ✅       |
   | 3.14   | ❌       | ✅       |

   = 6 cells × 4 runner OSes = **24**. The numpy1×3.13 and numpy1×3.14 cells are
   `exclude`d (no numpy 1.26 wheel). The 24 test cells reuse the 16 built wheels:
   the cp311/cp312 wheel for each OS is installed twice (once under numpy 1, once
   under numpy 2), which is what turns 16 wheels into 24 tested environments —
   no rebuild, just a second `pip install` with a different numpy constraint.

Boost is **not** needed in the test job — nothing rebuilds there. The extra YAML
over the fused job (artifact passing, picking the right wheel file per
Python/OS) is the price of not building each wheel twice.

## What "green" buys us

A fully green matrix is the signal to (a) drop `numpy < 2` in *AstroWISP*'s
`pyproject.toml` for real, then (b) drop `numpy < 2` and raise the
`requires-python` ceiling to allow 3.13/3.14 in *AutoWISP*'s `pyproject.toml`
(currently pinned `>=3.11,<3.13` precisely because astrowisp forces numpy 1).
