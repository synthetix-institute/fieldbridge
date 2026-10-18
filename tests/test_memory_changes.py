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


# ------------------------------------------------------------------------------------------------ carried over from
# the tests of the former `memory moves` (the same term handling, boundedness test, closure and isolated judgement)
import subprocess
import sys

PITCHFORK = json.loads((ROOT / "examples" / "memory" / "pitchfork.json").read_text(encoding="utf-8"))


def plain_pitchfork(drift="eps*x - x**3", control=("eps", [-0.5, 1.5]), **params):
    spec = dict(PITCHFORK)
    spec.pop("potential", None)
    spec["drift"] = {"x": drift}
    spec["parameters"] = dict(PITCHFORK["parameters"], **params)
    spec["control"] = {"name": control[0], "range": list(control[1])}
    spec["material"] = sorted(set(spec["parameters"]) - {control[0]})
    return spec


def test_terms_split_at_top_level_signs_only():
    assert ch.split_terms("eps*x - x**3 + h") == ["eps*x", "- x**3", "+ h"]
    assert ch.split_terms("-(K*(x**2 + y**2) - p + Delta)*x - 0.5*kappa*y") == [
        "-(K*(x**2 + y**2) - p + Delta)*x", "- 0.5*kappa*y"]
    assert ch.split_terms("2e-3*x + exp(-x)*y") == ["2e-3*x", "+ exp(-x)*y"]
    terms = ch.split_terms("a*x - b*x**3 + c")
    assert ch.join_terms(terms) == "a*x - b*x**3 + c" and ch.join_terms(terms[1:]) == "- b*x**3 + c"


def test_reasons_are_classed_from_the_derivation_texts():
    assert ch.reason_class("one-sided write (fold) and no material parameter reaches a cusp") == "one-sided write"
    assert ch.reason_class("positive cubic term: the write is symmetric but the state jumps") == "distant branch"
    assert ch.reason_class("degenerate critical mode: 2 directions become critical together") == "several critical directions"
    assert ch.reason_class("vanishing cubic term") == "vanishing cubic"
    assert ch.reason_class("no stable state (the realization oscillates)") == "oscillation"
    assert ch.reason_class("no oscillation: at 11 values of p every start settles") == "no oscillation"
    assert ch.reason_class("a single stable state: there is no second one") == "single stable state"
    assert ch.reason_class("something else entirely") == "unclassified"
    assert ch.reason_class("stopped at the time limit of 60 s") == "not derived"


def test_the_boundedness_test_keeps_the_pitchfork_and_rejects_it_without_its_cubic():
    ok, why = ch.bounded(plain_pitchfork())
    assert ok and why == ""
    ok, why = ch.bounded(plain_pitchfork("eps*x"))
    assert not ok and "unbounded" in why and "eps = 1.5" in why  # the growth is at the positive end of the range
    ok, why = ch.bounded(plain_pitchfork("eps*x + x**3"))        # the subcritical pitchfork without its saturation
    assert not ok


def test_the_four_kinds_of_change_and_the_donor_terms():
    spec = plain_pitchfork("eps*x - x**3 + h", h=0.1)
    rejected = []
    kept = {(k, w) for k, w, _ in ch.single_changes(spec, rejected=rejected)}
    assert ("control", "control = h over [0.025, 0.4]") in kept
    assert ("range", "range of eps wider: [-20.5, 21.5]") in kept
    assert ("remove", "remove '+ h' from dx/dt") in kept and ("remove", "remove 'eps*x' from dx/dt") in kept
    assert all(k != "add" for k, _ in kept)   # no donor
    dropped_cubic = [r for r in rejected if r["what"] == "remove '- x**3' from dx/dt"]
    assert dropped_cubic and "unbounded" in dropped_cubic[0]["why"]
    donor = {"name": "quintic donor", "drift": {"y": "a*y - b*y**5 + c*sin(y)"}, "parameters": {"a": 1.0, "b": 1.0, "c": 0.5}}
    added = [(w, s) for k, w, s in ch.single_changes(spec, donor=donor, check_bounded=False) if k == "add"]
    whats = {w for w, _ in added}
    assert "add g_add*(1*x) to dx/dt (from quintic donor)" in whats
    assert "add g_add*(1*x**5) to dx/dt (from quintic donor)" in whats
    assert "add g_add*(1*sin(x)) to dx/dt (from quintic donor)" in whats
    s5 = next(s for w, s in added if "x**5" in w)
    assert s5["parameters"]["g_add"] == 1.0 and "g_add" in s5["material"]


def test_a_combination_of_known_classes_closes():
    base = (("threshold-write",), (("symmetric-write", "one-sided write"),))
    other = ((), (("threshold-write", "single stable state"),))
    assert ch.closes(((), (("symmetric-write", "one-sided write"), ("threshold-write", "single stable state"))),
                     [base, other]) == (True, "a combination of known outcomes")


def test_the_isolated_judge_reads_the_pitchfork_and_the_fold_as_the_derivations_do():
    """Each derivation runs in a child process under a wall-clock limit; the result equals the in-process one."""
    j = ch.judge(plain_pitchfork(), budget=120)
    assert j["symmetric-write"]["class"] == "reached" and j["threshold-write"]["class"] == "reached"
    assert ch.judge(plain_pitchfork(), isolate=False)["symmetric-write"]["class"] == "reached"
    j = ch.judge(plain_pitchfork("eps*x + x**3 - x**5", control=("eps", [-0.5, 0.5])), budget=120)
    assert j["symmetric-write"]["status"] == "obstructed" and j["symmetric-write"]["class"] == "one-sided write"


def test_the_changes_command_with_a_donor(tmp_path):
    spec = plain_pitchfork("eps*x - x**3 + h", h=0.1)
    path = tmp_path / "biased.json"
    path.write_text(json.dumps(spec), encoding="utf-8")
    donor = tmp_path / "donor.json"
    # a donor term enters with a plus sign and the donor's parameters set to one; sin(x) keeps the states bounded
    donor.write_text(json.dumps({"name": "donor", "drift": {"y": "c*sin(y)"}, "parameters": {"c": 0.5}}), encoding="utf-8")
    out = tmp_path / "changes"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "changes", str(path), "--out-dir", str(out),
           "--no-isolate", "--kinds", "remove", "--donor", str(donor)]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    report = json.loads((out / "changes.json").read_text(encoding="utf-8"))
    assert report["command"] == "changes" and report["settings"]["donor"] == str(donor)
    kinds = {r["kind"] for r in report["singles"]}
    assert kinds == {"remove", "add"}
    assert any(r["what"] == "remove '- x**3' from dx/dt" for r in report["rejected"])
