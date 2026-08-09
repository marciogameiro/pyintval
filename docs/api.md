# API reference

All names below are available directly on the top-level `pyintval` module.

## Interval

```{eval-rst}
.. autoclass:: pyintval.Interval
   :members:
   :undoc-members:
   :special-members: __contains__, __and__, __or__
```

## DecoratedInterval

```{eval-rst}
.. autoclass:: pyintval.DecoratedInterval
   :members:
   :undoc-members:
```

## Elementary functions

Each math function is polymorphic: an `Interval` argument returns an `Interval`,
a `DecoratedInterval` returns a `DecoratedInterval` carrying the propagated
decoration.

```{eval-rst}
.. autofunction:: pyintval.sqrt
.. autofunction:: pyintval.sqr
.. autofunction:: pyintval.cbrt
.. autofunction:: pyintval.exp
.. autofunction:: pyintval.exp2
.. autofunction:: pyintval.exp10
.. autofunction:: pyintval.expm1
.. autofunction:: pyintval.log
.. autofunction:: pyintval.log2
.. autofunction:: pyintval.log10
.. autofunction:: pyintval.log1p
.. autofunction:: pyintval.sin
.. autofunction:: pyintval.cos
.. autofunction:: pyintval.tan
.. autofunction:: pyintval.asin
.. autofunction:: pyintval.acos
.. autofunction:: pyintval.atan
.. autofunction:: pyintval.atan2
.. autofunction:: pyintval.sinh
.. autofunction:: pyintval.cosh
.. autofunction:: pyintval.tanh
.. autofunction:: pyintval.asinh
.. autofunction:: pyintval.acosh
.. autofunction:: pyintval.atanh
.. autofunction:: pyintval.hypot
.. autofunction:: pyintval.pow
.. autofunction:: pyintval.pown
.. autofunction:: pyintval.erf
.. autofunction:: pyintval.erfc
```

## Rounding, sign, and step functions

```{eval-rst}
.. autofunction:: pyintval.abs
.. autofunction:: pyintval.fma
.. autofunction:: pyintval.recip
.. autofunction:: pyintval.floor
.. autofunction:: pyintval.ceil
.. autofunction:: pyintval.trunc
.. autofunction:: pyintval.round
.. autofunction:: pyintval.round_ties_to_away
.. autofunction:: pyintval.sign
.. autofunction:: pyintval.min
.. autofunction:: pyintval.max
```

## Set operations and reverse (constraint) operations

```{eval-rst}
.. autofunction:: pyintval.hull
.. autofunction:: pyintval.intersection
.. autofunction:: pyintval.cancel_minus
.. autofunction:: pyintval.cancel_plus
.. autofunction:: pyintval.mul_rev
.. autofunction:: pyintval.sqr_rev
.. autofunction:: pyintval.abs_rev
```

## Constructors and constants

```{eval-rst}
.. autofunction:: pyintval.empty
.. autofunction:: pyintval.entire
.. autofunction:: pyintval.pi
.. autofunction:: pyintval.e
```
