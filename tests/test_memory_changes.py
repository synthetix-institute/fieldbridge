"""Changes of a model and their composition (memory changes)."""
import json
from pathlib import Path

from fieldbridge.memory import changes as ch

ROOT = Path(__file__).resolve().parents[1]
HOPF = json.loads((ROOT / "examples/changes/hopf_normal_form.json").read_text(encoding="utf-8"))


def test_terms_split_and_join_back():
    assert ch.split_terms("-x**3 * omega - x * y**2 * omega + mu * omega * x - omega * y") == [
        "-x**3 * omega", "- x * y**2 * omega", "+ mu * omega * x", "- omega * y"]
    # exponents and brackets stay inside their term; a parameter named e is a term of its own
    assert ch.split_terms("1.5e-3*x - (a + b)*y") == ["1.5e-3*x", "- (a + b)*y"]
    assert ch.split_terms("a*x1 - x1**3 - e*(x2 + x3)") == ["a*x1", "- x1**3", "- e*(x2 + x3)"]
    assert ch.split_terms("x - e - y") == ["x", "- e", "- y"]
    assert ch.join_terms(["+ x", "- y", "z"]) == "x - y + z"


def test_a_change_that_removes_the_saturation_is_discarded():
    spec = json.loads((ROOT / "examples/memory/pitchfork.json").read_text(encoding="utf-8"))
    assert spec["drift"] == {"x": "eps*x - x**3"}
    rejected = []
    kept = ch.single_changes(spec, kinds=("remove",), rejected=rejected)
    assert [w for _, w, _ in kept] == ["remove 'eps*x' from dx/dt"]
    assert rejected and rejected[0]["what"] == "remove '- x**3' from dx/dt" and "unbounded" in rejected[0]["why"]


def test_closure_of_a_pair():
    base = ((), (("symmetric-write", "oscillation"), ("threshold-write", "single stable state")))
    one = ((), (("symmetric-write", "several critical directions"), ("threshold-write", "single stable state")))
    assert ch.closes(one, [base, one]) == (True, "equals a known outcome")
    ok, why = ch.closes((("symmetric-write",), (("threshold-write", "single stable state"),)), [base, one])
    assert not ok and "symmetric-write" in why
    ok, why = ch.closes(((), (("symmetric-write", "distant branch"), ("threshold-write", "single stable state"))),
                        [base, one])
    assert not ok and "distant branch" in why
    # a derivation that did not resolve decides nothing
    assert ch.closes(((), (("symmetric-write", "not derived"), ("threshold-write", "single stable state"))),
                     [base, one])[0]


def test_the_hopf_form_becomes_a_switch_only_when_two_terms_are_removed():
    """Removing the rotation term -omega y alone turns the rotation into a shear (two directions critical together);
    removing the growth term mu omega x alone leaves the oscillation; removing both leaves a pitchfork in y."""
    base = ch.judge(HOPF, budget=120)
    assert base["symmetric-write"]["class"] == "oscillation"
    singles = {w: s for _, w, s in ch.single_changes(HOPF, kinds=("remove",))}
    shear = singles["remove '- omega * y' from dx/dt"]
    growth = singles["remove '+ mu * omega * x' from dx/dt"]
    j1, j2 = ch.judge(shear, budget=120), ch.judge(growth, budget=120)
    assert j1["symmetric-write"]["class"] == "several critical directions"
    assert j2["symmetric-write"]["class"] == "oscillation"
    both = dict(shear)
    both["drift"] = dict(shear["drift"], x=ch.join_terms(
        [t for t in ch.split_terms(shear["drift"]["x"]) if t != "+ mu * omega * x"]))
    jb = ch.judge(both, budget=120)
    assert jb["symmetric-write"]["class"] == "reached"
    ok, why = ch.closes(ch.outcome(jb), [ch.outcome(base), ch.outcome(j1), ch.outcome(j2)])
    assert not ok and "symmetric-write" in why
