"""Milestone 1 smoke tests: the package imports and the extension is sound."""

import pyintval


def test_version_string():
    assert isinstance(pyintval.__version__, str)
    assert pyintval.__version__[0].isdigit()


def test_kernel_abi_version():
    assert pyintval.KERNEL_ABI_VERSION == 1


def test_build_info_reports_ieee754():
    info = pyintval.build_info()
    assert info["ieee754_doubles"] is True
    assert info["cxx_standard"] >= 202002
