"""Contributing a material: the template, the checks a contribution must pass, the catalog and the demo."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sympy")

from fieldbridge.memory import contribute, spec  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "memory"


@pytest.mark.parametrize("carrier", ["euclid", "orthant", "torus"])
def test_template_runs_but_is_not_ready_until_its_placeholders_are_replaced(carrier):
    template = contribute.template("my_material", carrier)
    spec.load(template)  # the placeholder dynamics loads and runs
    result = contribute.check(template, structure_only=True)
    failed = {r["check"] for r in result["checks"] if r["required"] and not r["passed"]}
    assert not result["passed"]
    assert failed == {"question", "source", "field", "assumptions"}
    reference = next(r for r in result["checks"] if r["check"] == "reference")
    assert not reference["passed"]  # the placeholder 'arXiv:0000.00000' is not a reference


def test_a_filled_template_is_ready_to_contribute():
    material = contribute.template("nematic_cell", "torus")
    material.update(question="Does a nematic cell with two easy axes store two orientations?",
                    field="soft matter", closure="thermal bath; anchoring fixed at the plates",
                    observable="birefringence between crossed polarizers",
                    assumptions=["one orientation angle per cell", "anchoring energy k in units of kT"])
    material["provenance"].update(source="test material written for this test (2026)",
                                  contributor="tester")
    result = contribute.check(material, structure_only=True)
    assert result["passed"]
    assert all(r["passed"] for r in result["checks"])


def test_every_example_material_is_ready_to_contribute():
    paths = contribute.specifications([EX])
    assert len(paths) >= 20
    for path in paths:
        result = contribute.check(path, structure_only=True)
        assert result["passed"], (path.name, [r for r in result["checks"] if r["required"] and not r["passed"]])


def test_check_reports_a_specification_that_does_not_load():
    broken = contribute.template("broken")
    broken["drift"] = {"x": "mu*y"}
    result = contribute.check(broken, structure_only=True)
    assert not result["passed"] and result["checks"][0]["check"] == "loads"
    assert "Undeclared symbol 'y'" in result["checks"][0]["message"]


def test_catalog_lists_each_material_with_its_calculated_memory():
    result = contribute.catalog([EX / "toggle.json", EX / "fields"], root=ROOT)
    assert [r["name"] for r in result["materials"]] == ["genetic toggle switch"]  # field specifications are skipped
    md = contribute.catalog_markdown(result)
    assert "[genetic toggle switch](../examples/memory/toggle.json)" in md
    assert "supercritical pitchfork" in md


def test_command_line_new_and_check(tmp_path):
    out = tmp_path / "cell.json"
    run = lambda *a: subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "memory", *a], cwd=ROOT,
                                    capture_output=True, text=True)
    assert run("new", "cell", "--carrier", "orthant", "--out", str(out)).returncode == 0
    assert json.loads(out.read_text(encoding="utf-8"))["carrier"]["kind"] == "orthant"
    assert run("new", "cell", "--out", str(out)).returncode != 0  # an existing file is not overwritten
    failing = run("check", str(out), "--structure-only")
    assert failing.returncode == 1 and "Not ready" in failing.stdout
    passing = run("check", str(EX / "toggle.json"), "--structure-only")
    assert passing.returncode == 0 and "Ready to contribute" in passing.stdout


def test_demo_writes_one_page_with_three_calculated_results(tmp_path):
    pytest.importorskip("matplotlib")
    from fieldbridge import demo
    results = demo.run(tmp_path, log=lambda *a: None)
    page = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert results["material"]["states"] == 2 and results["material"]["agreement"]
    assert abs(results["material"]["write_point"] - 2.0) < 0.05
    assert sum(str(r["status"]).startswith("reached") for r in results["mechanism"]["rows"]) == 7
    assert results["quantum"]["summary"]["reached"] == 6
    assert page.count("data:image/png;base64,") == 3 and "Adding a material specification" in page
    assert (tmp_path / "demo.md").exists()
