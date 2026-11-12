"""Accuracy regression: the quick self-test must pass (known geometry, exact answer)."""
from archaemetry.selftest import check, cases


def test_selftest_quick(tmp_path):
    res = cases(str(tmp_path))
    assert not check(res), check(res)
    assert all(r["verdict"] == "GO" for r in res)
    assert any(r.get("section_kind") == "section" for r in res)
