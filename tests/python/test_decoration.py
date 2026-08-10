"""DecoratedInterval: decoration propagation and the polymorphic math API."""

import pickle

import pytest

import pyintval as iv

DI = iv.DecoratedInterval


def test_construction_decorations():
    assert DI(1, 2).decoration == "com"
    assert DI(0, float("inf")).decoration == "dac"
    assert DI("[entire]").decoration == "dac"
    assert DI("[empty]").decoration == "trv"
    assert DI.nai().decoration == "ill"
    assert DI.nai().is_nai


def test_from_parts_and_properties():
    d = DI.from_parts(iv.Interval(1, 2), "def")
    assert d.decoration == "def"
    assert d.interval == iv.Interval(1, 2)
    assert d.lo == 1.0 and d.hi == 2.0
    assert d.is_defined and not d.is_defined_and_continuous and not d.is_common
    with pytest.raises(ValueError):
        DI.from_parts(iv.Interval(1, 2), "bogus")


def test_predicate_ladder():
    com = DI(1, 2)
    assert com.is_common and com.is_defined_and_continuous and com.is_defined
    dac = DI(0, float("inf"))
    assert not dac.is_common and dac.is_defined_and_continuous and dac.is_defined
    trv = iv.log(DI(0, 2))
    assert not trv.is_defined and not trv.is_defined_and_continuous


def test_continuous_ops_keep_com():
    a, b = DI(1, 2), DI(3, 4)
    assert (a + b).decoration == "com"
    assert (a * b).decoration == "com"
    assert iv.exp(a).decoration == "com"
    assert iv.sin(a).decoration == "com"
    assert iv.sqrt(DI(1, 4)).decoration == "com"


def test_domain_violations_give_trv():
    assert iv.sqrt(DI(-1, 4)).decoration == "trv"
    assert iv.log(DI(0, 2)).decoration == "trv"
    assert (DI(1, 2) / DI(-1, 1)).decoration == "trv"
    assert iv.asin(DI(-2, 0.5)).decoration == "trv"
    assert iv.tan(DI(1.5, 1.7)).decoration == "trv"
    assert iv.atan2(DI(-1, 1), DI(-1, 1)).decoration == "trv"  # origin


def test_step_functions_and_branch_cut():
    assert iv.floor(DI(0.2, 0.8)).decoration == "com"
    assert iv.floor(DI(0.5, 1.5)).decoration == "def"
    assert iv.sign(DI(1, 2)).decoration == "com"
    assert iv.sign(DI(-1, 1)).decoration == "def"
    assert iv.atan2(DI(-1, 1), DI(-2, -1)).decoration == "def"  # branch cut


def test_composition_propagates_weakest():
    # A single domain violation poisons the whole certificate.
    assert iv.exp(iv.sqrt(DI(-1, 4))).decoration == "trv"
    assert (iv.log(DI(0, 2)) + DI(1, 2)).decoration == "trv"
    # A clean composition certifies defined & continuous.
    f = iv.sin(DI(0, 1)) / (iv.sqr(DI(0, 1)) + DI(1, 1))
    assert f.is_defined_and_continuous
    # ill poisons regardless.
    assert (DI.nai() + DI(1, 2)).decoration == "ill"


def test_operators_and_promotion():
    a = DI(2, 3)
    assert (a + 1).decoration == "com"
    assert (2 * a).interval == iv.Interval(4, 6)
    assert (a + iv.Interval(1, 1)).decoration == "com"  # bare interval promoted
    assert (a**2).interval == iv.Interval(4, 9)
    assert 2**0.5 in (DI(2) ** 0.5).interval
    assert (-a).interval == iv.Interval(-3, -2)


def test_polymorphism_preserves_type():
    assert isinstance(iv.exp(iv.Interval(0, 1)), iv.Interval)
    assert isinstance(iv.exp(DI(0, 1)), DI)
    assert isinstance(iv.atan2(DI(1, 2), DI(1, 2)), DI)
    assert isinstance(iv.hypot(iv.Interval(3), iv.Interval(4)), iv.Interval)


def test_repr_str_pickle():
    a = DI(1, 2)
    assert "com" in repr(a)
    assert str(a).endswith("_com")
    assert pickle.loads(pickle.dumps(a)) == a
    assert eval(repr(a), {"DecoratedInterval": DI, "Interval": iv.Interval}) == a


def test_decorated_set_cancel_reverse_ops():
    # F3: set operations, cancellative subtraction, and reverse ops accept
    # decorated intervals; the result is the bare op decorated trv, and NaI
    # propagates.
    I = iv.Interval
    fp = DI.from_parts
    a, b = fp(I(1, 3), "com"), fp(I(2.1, 4), "com")
    assert iv.intersection(a, b) == fp(I(2.1, 3), "trv")
    assert iv.hull(a, b) == fp(I(1, 4), "trv")
    assert iv.cancel_plus(fp(I(-1, 3), "com"), fp(I(0, 1), "com")).decoration == "trv"
    assert iv.sqr_rev(fp(I(0, 25), "def")) == fp(I(-5, 5), "trv")
    assert iv.mul_rev(a, b).decoration == "trv"
    assert iv.abs_rev(fp(I(-1, 2), "com")).decoration == "trv"
    # NaI poisons the result.
    assert iv.intersection(DI.nai(), a).is_nai
    assert iv.mul_rev(DI.nai(), b).is_nai
    # Bare Interval overloads still work unchanged.
    assert iv.intersection(I(1, 3), I(2, 4)) == I(2, 3)


def test_decorated_text_literals():
    # F4: native "[a,b]_dec" parsing; no suffix -> the tightest decoration; the
    # uncertain form works decorated too.
    assert DI("[1,1e3]_com") == DI.from_parts(iv.Interval(1, 1000), "com")
    assert DI("[1,1E3]_COM").decoration == "com"  # case-insensitive suffix
    assert DI("[entire]").decoration == "dac"
    assert DI("3.56?1e2") == DI.from_parts(iv.Interval(355, 357), "com")
    # Invalid decoration -> NaI (over-claim, unknown, or 'ill'); "[nai]" -> NaI.
    assert DI("[nai]").is_nai
    assert DI("[empty]_com").is_nai
    assert DI("[,]_com").is_nai
    assert DI("[1,2]_ill").is_nai
    assert DI("[1,2]_fooo").is_nai
    # Genuinely malformed still raises.
    with pytest.raises(ValueError):
        DI("garbage")
