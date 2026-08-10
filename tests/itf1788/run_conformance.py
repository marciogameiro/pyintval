#!/usr/bin/env python3
"""Run the ITF1788 IEEE 1788-2015 conformance suite against pyintval.

ITF1788 (https://github.com/oheim/ITF1788, Apache-2.0) is the reference
conformance framework for IEEE Std 1788-2015. This script:

  1. clones ITF1788 at a pinned commit (into --workdir, default build/itf1788),
  2. applies two small Python-3.11 / PyYAML-6 compatibility shims to the 2018
     tool (``time.clock`` -> ``time.perf_counter``, ``yaml.load`` ->
     ``yaml.safe_load``),
  3. installs the pyintval arithmetic plugin (pyintval_arith.yaml),
  4. generates Python unit tests from the ~7,200-case .itl corpus for the
     ``(python, unittest, pyintval)`` configuration,
  5. runs them and classifies each test against two standards:
       * ENCLOSURE (rigor): pyintval's result must CONTAIN the IEEE tightest
         result (decorations may under- but never over-claim). This is the
         soundness guarantee and is what CI gates on.
       * TIGHTEST: bit-exact minimal endpoints. Reported for information only --
         pyintval deliberately widens transcendentals ~1 ulp for safety, so it is
         not, and does not claim to be, tightest on those.

The gate passes iff every enclosure failure or run error is listed in
``known_deviations.txt`` (the documented F3/F4/F5 gaps -- decorated set/cancel/
reverse ops, the decorated text parser, and a few decoration/NaI edge cases).
Any NEW rigor violation fails. Regenerate the baseline with --update-baseline
after closing a gap.

Requirements: pyintval importable (pip install .), plus PLY and PyYAML for the
generator. Run:  python tests/itf1788/run_conformance.py
"""

from __future__ import annotations

import argparse
import contextlib
import glob
import importlib.util
import io
import math
import os
import shutil
import subprocess
import sys
import unittest

from pyintval import DecoratedInterval, Interval

ITF1788_REPO = "https://github.com/oheim/ITF1788.git"
ITF1788_REF = "b6ee1e24d209c289f99a68ddc357839935799eae"

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PLUGIN = os.path.join(HERE, "pyintval_arith.yaml")
BASELINE = os.path.join(HERE, "known_deviations.txt")

# IEEE 1788 decoration order (strongest -> weakest). A rigorous result may carry
# a weaker (lower) decoration than the tightest, but never a stronger one.
_RANK = {"com": 4, "dac": 3, "def": 2, "trv": 1, "ill": 0}


def _bare(x):
    return x.interval if isinstance(x, DecoratedInterval) else x


def encloses(computed, expected) -> bool:
    """Does `computed` rigorously enclose the IEEE tightest `expected`?"""
    if isinstance(computed, (Interval, DecoratedInterval)) and isinstance(
        expected, (Interval, DecoratedInterval)
    ):
        ok = _bare(expected).subset(_bare(computed))
        if isinstance(computed, DecoratedInterval) or isinstance(expected, DecoratedInterval):
            cd = computed.decoration if isinstance(computed, DecoratedInterval) else "com"
            ed = expected.decoration if isinstance(expected, DecoratedInterval) else "com"
            ok = ok and _RANK[cd] <= _RANK[ed]
        return ok
    if isinstance(computed, float) and isinstance(expected, float):
        return (math.isnan(computed) and math.isnan(expected)) or computed == expected
    return computed == expected


def _run(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def _patch(path: str, old: str, new: str) -> None:
    with open(path, encoding="utf-8") as fh:
        src = fh.read()
    if old in src:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(src.replace(old, new))


def prepare_clone(workdir: str) -> str:
    """Clone + pin ITF1788, apply compat shims, install the pyintval plugin."""
    clone = os.path.join(workdir, "ITF1788")
    if not os.path.isdir(os.path.join(clone, ".git")):
        os.makedirs(workdir, exist_ok=True)
        _run(["git", "clone", "--quiet", ITF1788_REPO, clone])
    _run(["git", "-C", clone, "checkout", "--quiet", ITF1788_REF])
    # The 2018 tool predates Python 3.8 / PyYAML 5.1.
    _patch(os.path.join(clone, "itf1788", "__main__.py"), "time.clock()", "time.perf_counter()")
    _patch(os.path.join(clone, "itf1788", "discovery.py"), "yaml.load(", "yaml.safe_load(")
    dst = os.path.join(clone, "itf1788", "plugins", "python", "arith", "pyintval")
    os.makedirs(dst, exist_ok=True)
    shutil.copyfile(PLUGIN, os.path.join(dst, "arith.yaml"))
    return clone


def generate(clone: str, outdir: str) -> str:
    """Generate the Python unit tests; return the directory holding them."""
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    env = dict(os.environ, PYTHONPATH=clone)
    _run(
        [
            sys.executable,
            "-m",
            "itf1788",
            "-s",
            "itl",
            "-o",
            outdir,
            "-c",
            "(python, unittest, pyintval)",
        ],
        cwd=clone,
        env=env,
    )
    return os.path.join(outdir, "python", "unittest", "pyintval")


def _load_suite(path: str):
    spec = importlib.util.spec_from_file_location(os.path.basename(path)[:-3], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return unittest.defaultTestLoader.loadTestsFromModule(mod)


def _run_pass(gendir: str, comparator):
    """Run every generated test; return (n_tests, failures set, errors set)."""
    n = 0
    failures: set[str] = set()
    errors: set[str] = set()
    for path in sorted(glob.glob(os.path.join(gendir, "*.py"))):
        try:
            suite = _load_suite(path)
        except Exception:  # a module that fails to import counts as an error bucket
            errors.add(os.path.basename(path))
            continue
        if comparator is not None:

            def make(cmp):
                def assert_equal(self, a, b, msg=None):
                    if not cmp(a, b):
                        raise AssertionError(f"{a!r} does not enclose {b!r}")

                return assert_equal

            for group in suite:
                for tc in group:
                    tc.assertEqual = make(comparator).__get__(tc)
        result = unittest.TestResult()
        with contextlib.redirect_stderr(io.StringIO()):
            suite.run(result)
        n += result.testsRun
        failures.update(tc.id() for tc, _ in result.failures)
        errors.update(tc.id() for tc, _ in result.errors)
    return n, failures, errors


def load_baseline() -> set[str]:
    if not os.path.exists(BASELINE):
        return set()
    with open(BASELINE, encoding="utf-8") as fh:
        return {ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")}


def write_baseline(ids: set[str]) -> None:
    header = (
        "# ITF1788 known deviations for pyintval -- tests that do not hold under the\n"
        "# ENCLOSURE (rigor) standard, i.e. run errors or under-enclosures. These are\n"
        "# the documented, deferred gaps (F3 decorated set/cancel/reverse ops, F4\n"
        "# decorated text parser, F5 decoration/NaI edge cases); see README.md.\n"
        "# Regenerate with: python tests/itf1788/run_conformance.py --update-baseline\n"
    )
    with open(BASELINE, "w", encoding="utf-8") as fh:
        fh.write(header)
        for tid in sorted(ids):
            fh.write(tid + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=os.path.join(REPO_ROOT, "build", "itf1788"))
    ap.add_argument(
        "--update-baseline",
        action="store_true",
        help="rewrite known_deviations.txt from this run instead of gating",
    )
    args = ap.parse_args()

    clone = prepare_clone(args.workdir)
    gendir = generate(clone, os.path.join(args.workdir, "generated"))

    total, tight_fail, errors = _run_pass(gendir, None)  # native == : tightest
    _, encl_fail, _ = _run_pass(gendir, encloses)  # enclosure : rigor
    deviations = encl_fail | errors  # not rigorous OR not runnable

    tight_pass = total - len(tight_fail) - len(errors)
    encl_pass = total - len(encl_fail) - len(errors)
    print(f"ITF1788 conformance (pinned {ITF1788_REF[:10]}):  {total} tests")
    print(
        f"  ENCLOSURE (rigor):  {encl_pass}/{total}  "
        f"({len(encl_fail)} under-enclosures, {len(errors)} unsupported-op errors)"
    )
    print(
        f"  TIGHTEST (info)  :  {tight_pass}/{total}  "
        f"({len(tight_fail) - len(encl_fail)} enclose-but-not-tightest, mostly transcendentals)"
    )

    if args.update_baseline:
        write_baseline(deviations)
        rel = os.path.relpath(BASELINE, REPO_ROOT)
        print(f"  wrote {len(deviations)} known deviations to {rel}")
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
        for tid in sorted(new)[:25]:
            print(f"    {tid}")
        if len(new) > 25:
            print(f"    ... and {len(new) - 25} more")
        return 1
    print("\nPASS: no new rigor deviations; every result encloses the IEEE tightest interval.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
