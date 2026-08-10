"""Construction, parsing, error handling, repr/str, and pickling."""

import copy
import math
import pickle

import pytest

import pyintval as iv


def test_from_bounds():
    x = iv.Interval(1.0, 2.0)
    assert x.lo == 1.0
    assert x.hi == 2.0


def test_point_from_number():
    for v in (0, 1, -3, 2.5, 1e300):
        x = iv.Interval(v)
        assert x.lo == float(v)
        assert x.hi == float(v)
        assert x.is_singleton


def test_from_string_is_outward_rounded():
    x = iv.Interval("0.1")
    assert x.lo < x.hi  # 0.1 is not representable, so the enclosure is nondegenerate
    assert x.lo < 0.1 or x.lo == 0.1
    assert 0.1 <= x.hi
    # exactly-representable literal stays a point
    assert iv.Interval("0.5").is_singleton


def test_from_string_bracket_forms():
    assert iv.Interval("[1, 2]").endpoints == (1.0, 2.0)
    assert iv.Interval("[entire]").is_entire
    assert iv.Interval("[empty]").is_empty
    assert iv.Interval("[]").is_empty
    assert iv.Interval("[3, inf]").hi == math.inf


def test_from_string_uncertain_form():
    # IEEE 1788 "d?ruE" uncertain notation (F4).
    assert iv.Interval("-10?").endpoints == (-10.5, -9.5)  # half-ulp default radius
    assert iv.Interval("-10?u").endpoints == (-10.0, -9.5)  # one-sided up
    assert iv.Interval("-10?12").endpoints == (-22.0, 2.0)  # radius of 12 ulps
    assert iv.Interval("3.56?1e2").endpoints == (355.0, 357.0)  # trailing exponent
    assert iv.Interval("0.0?u").lo == 0.0
    x = iv.Interval("3.56?1")  # [3.55, 3.57], outward-rounded
    assert x.lo <= 3.55 and 3.57 <= x.hi


def test_from_string_half_bounded():
    assert iv.Interval("[,]").is_entire
    assert iv.Interval("[1,]").endpoints == (1.0, math.inf)
    assert iv.Interval("[,2]").endpoints == (-math.inf, 2.0)


def test_empty_and_entire_constructors():
    assert iv.Interval.empty().is_empty
    assert iv.Interval.entire().is_entire
    assert iv.empty().is_empty
    assert iv.entire().is_entire


@pytest.mark.parametrize("lo,hi", [(2.0, 1.0), (math.nan, 1.0), (math.inf, math.inf)])
def test_invalid_bounds_raise(lo, hi):
    with pytest.raises(ValueError):
        iv.Interval(lo, hi)


@pytest.mark.parametrize("s", ["nonsense", "[2,1]", "1.2.3", "[inf,inf]", "nan", ""])
def test_malformed_string_raises(s):
    with pytest.raises(ValueError):
        iv.Interval(s)


def test_infinite_point_raises():
    with pytest.raises(ValueError):
        iv.Interval(math.inf)


def test_repr_round_trips_exactly():
    for x in [
        iv.Interval("0.1"),
        iv.Interval(1, 2),
        iv.Interval.empty(),
        iv.Interval.entire(),
        iv.Interval(-3.5, 3.5),
        iv.Interval(2, math.inf),
    ]:
        assert eval(repr(x), {"Interval": iv.Interval}) == x


def test_str_is_readable_and_encloses():
    x = iv.Interval("0.1")
    s = str(x)
    assert s.startswith("[") and s.endswith("]")
    reparsed = iv.Interval(s)
    assert x.subset(reparsed)  # I/O rounds outward


@pytest.mark.parametrize(
    "x",
    [
        iv.Interval(1, 2),
        iv.Interval("0.1"),
        iv.Interval.empty(),
        iv.Interval.entire(),
        iv.Interval(-math.inf, 5.0),
    ],
)
def test_pickle_and_copy(x):
    assert pickle.loads(pickle.dumps(x)) == x
    assert copy.copy(x) == x
    assert copy.deepcopy(x) == x


def test_hash_consistency():
    a = iv.Interval(1, 2)
    b = iv.Interval(1.0, 2.0)
    assert a == b
    assert hash(a) == hash(b)
    d = {a: "x"}
    assert d[b] == "x"
