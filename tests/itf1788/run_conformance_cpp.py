#!/usr/bin/env python3
"""Run the ITF1788 IEEE 1788-2015 conformance suite against pyintval's C++ kernel.

This is the C++ counterpart of ``run_conformance.py``. Where that script drives
the corpus through the Python extension, this one compiles the generated tests
against the header-only C++ API directly (``include/pyintval/*.hpp`` + the
vendored CORE-MATH kernels), with no Python or pybind11 in the loop -- the same
conformance bar applied to the layer that CMGDB and other C++ callers reuse.

Pipeline:

  1. clone + pin ITF1788 and apply the Py-3.11 / PyYAML-6 shims (shared with the
     Python runner via ``run_conformance.prepare_clone``),
  2. install the C++ plugins -- ``pyintval_cpp_arith.yaml`` (op -> C++ API map)
     and ``cpp_doctest_test.yaml`` (render with doctest, already vendored), and
     stage the core corpus (see ``CORE_ITL``) into ``itl_core/``,
  3. generate the ``(cpp, doctest, pyintval)`` tests,
  4. build the CORE-MATH kernels once into a static archive,
  5. compile + link + run each generated file, gating on ENCLOSURE: every
     assertion asks whether pyintval's result CONTAINS the IEEE tightest result
     (``itf_eq`` in the arith plugin). A failing assertion is an under-enclosure.

The gate passes iff no new under-enclosure appears versus ``known_deviations_cpp.txt``
(empty: the C++ kernel encloses on the entire core corpus). Regenerate with
--update-baseline.

Excluded from the C++ target: the four corpus files exercising operations
pyintval does not implement (Allen ``overlap``; ``dot``/``sum`` reductions; the
reverse trigonometric ops ``sinRev``/``cosRev``/...; ``powRev``) and the three
vendor-extension files (c-xsc, fi_lib, mpfi) that use reciprocal trig
(``sec``/``csc``/``cot``/...) and ``rootn``. Those operations lie outside the
IEEE 1788 set pyintval provides; the Python runner covers the same core corpus.

Requirements: a C++20 compiler and a C compiler (``$CXX``/``$CC``, default
``c++``/``cc``), plus PLY and PyYAML for the generator. pyintval need NOT be
installed -- this target compiles the headers, not the wheel.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_conformance as rc  # noqa: E402  (shared clone/shim scaffolding + pinned ref)

REPO_ROOT = rc.REPO_ROOT
ARITH = os.path.join(HERE, "pyintval_cpp_arith.yaml")
TESTLIB = os.path.join(HERE, "cpp_doctest_test.yaml")
BASELINE = os.path.join(HERE, "known_deviations_cpp.txt")

# The core IEEE 1788 corpus files, exercising only operations pyintval implements.
CORE_ITL = [
    "abs_rev",
    "atan2",
    "ieee1788-constructors",
    "ieee1788-exceptions",
    "libieeep1788_bool",
    "libieeep1788_cancel",
    "libieeep1788_class",
    "libieeep1788_elem",
    "libieeep1788_mul_rev",
    "libieeep1788_num",
    "libieeep1788_rec_bool",
    "libieeep1788_set",
]

_ERR = re.compile(r"([^/\s]+\.cpp):(\d+): ERROR:")
_CASES = re.compile(r"test cases:\s*(\d+)")
_ASSERTS = re.compile(r"assertions:\s*(\d+)\s*\|\s*(\d+) passed\s*\|\s*(\d+) failed")


def prepare_clone(workdir: str) -> str:
    """Clone + shim ITF1788, install the C++ plugins, stage the core corpus."""
    # clone, pin, compat shims (also drops the Python plugin, harmless here):
    clone = rc.prepare_clone(workdir)
    tdst = os.path.join(clone, "itf1788", "plugins", "cpp", "test", "doctest")
    adst = os.path.join(clone, "itf1788", "plugins", "cpp", "arith", "pyintval")
    os.makedirs(tdst, exist_ok=True)
    os.makedirs(adst, exist_ok=True)
    shutil.copyfile(TESTLIB, os.path.join(tdst, "test.yaml"))
    shutil.copyfile(ARITH, os.path.join(adst, "arith.yaml"))
    coredir = os.path.join(clone, "itl_core")
    if os.path.isdir(coredir):
        shutil.rmtree(coredir)
    os.makedirs(coredir)
    for name in CORE_ITL:
        shutil.copyfile(
            os.path.join(clone, "itl", name + ".itl"), os.path.join(coredir, name + ".itl")
        )
    return clone


def generate(clone: str, outdir: str) -> str:
    """Generate the C++ conformance tests; return the directory holding them."""
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    env = dict(os.environ, PYTHONPATH=clone)
    rc._run(
        [
            sys.executable,
            "-m",
            "itf1788",
            "-s",
            "itl_core",
            "-o",
            outdir,
            "-c",
            "(cpp, doctest, pyintval)",
        ],
        cwd=clone,
        env=env,
    )
    return os.path.join(outdir, "cpp", "doctest", "pyintval")


def build_coremath(builddir: str, cc: str) -> str:
    """Compile the vendored CORE-MATH kernels once into a static archive."""
    os.makedirs(builddir, exist_ok=True)
    objs = []
    for src in sorted(glob.glob(os.path.join(REPO_ROOT, "third_party", "core-math", "*", "*.c"))):
        stem = os.path.splitext(os.path.basename(src))[0]
        tag = os.path.basename(os.path.dirname(src)) + "_" + stem
        obj = os.path.join(builddir, tag + ".o")
        rc._run([cc, "-O1", "-c", src, "-o", obj])
        objs.append(obj)
    lib = os.path.join(builddir, "libcoremath.a")
    rc._run(["ar", "rcs", lib, *objs])
    return lib


def compile_and_run(cpp: str, cxx: str, lib: str, builddir: str) -> dict:
    """Compile+link+run one generated file; return its counts and any deviations."""
    name = os.path.splitext(os.path.basename(cpp))[0]
    exe = os.path.join(builddir, "t_" + name)
    flags = [
        "-std=c++20",
        "-O0",
        "-I",
        os.path.join(REPO_ROOT, "include"),
        "-I",
        os.path.join(REPO_ROOT, "third_party"),
    ]
    cc = subprocess.run([cxx, *flags, cpp, lib, "-lm", "-o", exe], capture_output=True, text=True)
    if cc.returncode != 0:
        return {
            "name": name,
            "cases": 0,
            "assertions": 0,
            "failed": 1,
            "devs": {f"{name}.cpp:compile"},
            "log": (cc.stdout + cc.stderr),
        }
    rp = subprocess.run([exe, "--no-intro", "--no-version"], capture_output=True, text=True)
    out = rp.stdout + rp.stderr
    cases = int(m.group(1)) if (m := _CASES.search(out)) else 0
    asserts, failed = (int(a.group(1)), int(a.group(3))) if (a := _ASSERTS.search(out)) else (0, 1)
    return {
        "name": name,
        "cases": cases,
        "assertions": asserts,
        "failed": failed,
        "devs": {f"{f}:{ln}" for f, ln in _ERR.findall(out)},
        "log": out,
    }


def load_baseline() -> set[str]:
    if not os.path.exists(BASELINE):
        return set()
    with open(BASELINE, encoding="utf-8") as fh:
        return {ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")}


def write_baseline(ids: set[str]) -> None:
    header = (
        "# ITF1788 C++ target known deviations for pyintval -- generated tests that\n"
        "# do not hold under the ENCLOSURE (rigor) standard (a failing doctest CHECK,\n"
        "# keyed <file>:<line>). Empty: the C++ kernel encloses on the entire core\n"
        "# corpus. Regenerate with:\n"
        "#   python tests/itf1788/run_conformance_cpp.py --update-baseline\n"
    )
    with open(BASELINE, "w", encoding="utf-8") as fh:
        fh.write(header)
        for tid in sorted(ids):
            fh.write(tid + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=os.path.join(REPO_ROOT, "build", "itf1788_cpp"))
    ap.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    ap.add_argument("--cc", default=os.environ.get("CC", "cc"))
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument(
        "--update-baseline",
        action="store_true",
        help="rewrite known_deviations_cpp.txt from this run instead of gating",
    )
    args = ap.parse_args()

    clone = prepare_clone(args.workdir)
    gendir = generate(clone, os.path.join(args.workdir, "generated"))
    builddir = os.path.join(args.workdir, "build")
    if os.path.isdir(builddir):
        shutil.rmtree(builddir)
    lib = build_coremath(builddir, args.cc)

    files = sorted(glob.glob(os.path.join(gendir, "*.cpp")))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda f: compile_and_run(f, args.cxx, lib, builddir), files))

    total_cases = sum(r["cases"] for r in results)
    total_asrt = sum(r["assertions"] for r in results)
    deviations: set[str] = set()
    for r in sorted(results, key=lambda r: r["name"]):
        deviations |= r["devs"]
    print(
        f"ITF1788 C++ conformance (pinned {rc.ITF1788_REF[:10]}):  "
        f"{total_cases} test cases, {total_asrt} assertions across {len(files)} core files"
    )
    print(
        f"  ENCLOSURE (rigor):  {total_asrt - len(deviations)}/{total_asrt}  "
        f"({len(deviations)} under-enclosures / compile failures)"
    )
    print("  Corpus: 12 core IEEE 1788 files; 7 excluded (ops outside pyintval's set; see README).")

    if args.update_baseline:
        write_baseline(deviations)
        print(
            f"  wrote {len(deviations)} known deviations to "
            f"{os.path.relpath(BASELINE, REPO_ROOT)}"
        )
        return 0

    baseline = load_baseline()
    new = deviations - baseline
    stale = baseline - deviations
    if stale:
        print(
            f"  note: {len(stale)} baselined deviations now pass -- "
            f"run --update-baseline to prune them."
        )
    if new:
        print(f"\nFAIL: {len(new)} NEW rigor deviation(s) not in the baseline:")
        for r in results:
            offenders = r["devs"] & new
            if offenders:
                print(f"    {r['name']}: {', '.join(sorted(offenders))}")
                print("\n".join("      " + ln for ln in r["log"].splitlines()[:6]))
        return 1
    print(
        "\nPASS: no new rigor deviations; the C++ kernel encloses the IEEE tightest "
        "interval on the entire core corpus."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
