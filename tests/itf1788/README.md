# ITF1788 conformance suite

[ITF1788](https://github.com/oheim/ITF1788) is the reference conformance test
framework for **IEEE Std 1788-2015** (interval arithmetic). It compiles a
~7,200-case `.itl` test corpus — whose expected results are the *tightest*
correctly-rounded binary64 intervals, written as exact hex-float endpoints — into
unit tests for a target library. This directory wires pyintval into it.

## Running

```bash
pip install .                              # pyintval must be importable
pip install PLY PyYAML                      # the generator's dependencies
python tests/itf1788/run_conformance.py     # clone, generate, run, gate
```

The runner clones ITF1788 at a pinned commit into `build/itf1788/`, generates the
`(python, unittest, pyintval)` tests, runs them, and gates. It needs network
access on the first run only (to clone). Exit code is 0 on success.

## What it checks — two standards

Every generated test compares pyintval's result to the IEEE *tightest* interval.
We score each against two bars:

- **Enclosure (rigor)** — pyintval's result must **contain** the tightest
  interval (and a decorated result may carry a weaker decoration but never a
  stronger one). This is the soundness guarantee. **CI gates on this.**
- **Tightest** — bit-exact minimal endpoints. Reported for information only.

Current status (pinned corpus, 7,236 tests):

| Standard | Pass | |
|---|---|---|
| **Enclosure (rigor)** | **7,056 / 7,236 (98%)** | 99.5% of the tests that exercise an *implemented* operation |
| Tightest | 5,570 / 7,236 (77%) | the ~1,500 gap is transcendentals, ~1 ulp wider than tightest **by design** |

pyintval builds its transcendentals on CORE-MATH's correctly-rounded kernels and
**deliberately widens them one ulp per endpoint** for safety, so it is not — and
does not claim to be — *tightest* on `exp/log/trig/pow/…`. It remains a rigorous
enclosure, which is what matters and what the gate enforces.

## The gate and the baseline

`run_conformance.py` fails only on a **new** rigor deviation — an under-enclosure
or a run error not already recorded in [`known_deviations.txt`](known_deviations.txt).
That baseline (180 entries) captures the documented, deferred gaps below. After
closing one, prune the baseline:

```bash
python tests/itf1788/run_conformance.py --update-baseline
```

## Known gaps (deferred, tracked as follow-ups)

The suite surfaced four real gaps beyond the by-design transcendental widening
(F1); **F2 and F3 are now fixed**, F4 and F5 remain (baselined). None is a
soundness bug in an implemented operation.

| ID | Gap | ~Count | Baseline bucket |
|---|---|---|---|
| **F4** | `DecoratedInterval("[a,b]_com")` and some numeric/text literal forms aren't parsed by pyintval's decorated text constructor. | ~130 | `text_to_*`, `nums_to_*`, `IEEE1788_c/e` |
| **F5** | Decoration edge cases on discontinuous functions (`ceil/floor/trunc/round`) and NaI comparison semantics. | ~35 | `*_dec_test` under-enclosures |

**F2 (fixed):** `pow(x, y)` with a base touching 0 previously returned `[empty]`
via `exp(y·log 0)`; it now implements the IEEE 1788 corner semantics
(`elementary.hpp`), so those ~180 cases enclose.

**F3 (fixed):** `intersection`, `hull`, `cancel_minus/plus`, and the reverse ops
`mul_rev/sqr_rev/abs_rev` now have `DecoratedInterval` overloads (result `trv`,
NaI propagating; `decoration.hpp` + bindings), clearing ~180 run errors.
Regression tests for both live in `tests/cpp/test_{elementary,decoration}.cpp`
and `tests/python/test_{regressions,decoration}.py`.

## How the plugin works

`pyintval_arith.yaml` is the ITF1788 *arithmetic-library plugin*: it maps each
IEEE 1788 operation name to pyintval's Python API. The language/test-framework
plugins (Python + `unittest`) are ITF1788's own. Notes:

- Numeric functions and comparisons unwrap a decorated operand (`_iv`), since
  IEEE 1788 defines them on the interval part.
- A `_dtext` adapter parses the `[a,b]_dec` decorated-literal syntax that
  pyintval's own string constructor does not accept (this is F4; the adapter
  lets the *semantic* tests run anyway).
- Exception predicates are left empty: pyintval is set-based and never signals,
  so ITF1788's `signal` clauses become no-ops and only the interval result is
  checked.

## Provenance

ITF1788 is Apache-2.0 (© 2014 Nehmeier & Kiesner; plugins © 2015 Heimlich). It is
pinned at commit `b6ee1e2` and used as an unmodified dev-time tool — only cloned,
never vendored — with two in-place compatibility shims applied by the runner
(`time.clock` → `time.perf_counter`, `yaml.load` → `yaml.safe_load`) so the 2018
codebase runs on Python 3.11 / PyYAML 6. Only `pyintval_arith.yaml`,
`run_conformance.py`, and `known_deviations.txt` live in this repo.
