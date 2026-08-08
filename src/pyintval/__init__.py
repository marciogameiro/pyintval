"""Rigorous interval arithmetic for Python.

pyintval implements IEEE 1788-2015 set-based interval arithmetic with
correctly rounded double-precision endpoints: every operation returns an
interval mathematically guaranteed to contain the true result.

Milestone 1: packaging skeleton. The interval types and arithmetic land in
subsequent milestones.
"""

from __future__ import annotations

from pyintval._core import KERNEL_ABI_VERSION, __version__, build_info

__all__ = ["KERNEL_ABI_VERSION", "__version__", "build_info"]
