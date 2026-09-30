"""The instrument of the web page, run in Node against a stubbed document: changes of one component move along the
verified edges, a sequence applies its steps in order, a realization is reached by the shortest path of changes, and
the consequences shown are those calculated for the current realization."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sympy")

from fieldbridge import site_data as sd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NODE = shutil.which("node")
pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")
SUBSET = ["pitchfork", "pitchfork_below", "pitchfork_bias", "pitchfork_subcritical", "two_spins", "two_spins_x1",
          "two_spins_both", "two_spins_h0", "two_spins_z0", "spin1_transverse", "spin1_easy_axis", "stoner_wohlfarth", "sw_oblique", "sw_easy"]


@pytest.fixture(scope="module")
def data_file(tmp_path_factory):
    data = sd.build(only=SUBSET, derive=False, strict=False, log=lambda *a: None)  # the navigation needs no derivation
    return sd.write(data, tmp_path_factory.mktemp("site") / "data.js")


def scenario(data_file, steps, tmp_path):
    f = tmp_path / "scenario.json"
    f.write_text(json.dumps(steps), encoding="utf-8")
    done = subprocess.run([NODE, str(ROOT / "tests" / "site_ui_harness.cjs"), str(data_file), str(f)],
                          capture_output=True, text=True, check=True)
    return json.loads(done.stdout)


def test_the_page_opens_on_a_mechanism_written_without_a_field(data_file, tmp_path):
    first = scenario(data_file, [{"do": "look"}], tmp_path)[0]
    assert first["node"] == "pitchfork" and first["path"] == ["pitchfork"]
    assert first["mechanism"] == "symmetric write (supercritical pitchfork)" and first["lead"] == "2 stable states"
    # every change names the mechanism it leads to
    for name in ("one state", "one-sided write", "distant write", "same mechanism in magnetism"):
        assert name in first["changes"], name


def test_one_change_of_the_normal_form_gives_another_mechanism(data_file, tmp_path):
    got = scenario(data_file, [{"do": "edge", "edge": "pf_below"}, {"do": "undo"}, {"do": "edge", "edge": "pf_bias"},
                               {"do": "undo"}, {"do": "edge", "edge": "pf_subcritical"}], tmp_path)
    assert (got[0]["mechanism"], got[0]["lead"]) == ("single stable state", "one stable state")
    # the selection is shown beside the expression at the top of the page
    assert got[0]["selected"] == "single stable state | pitchfork normal form below the transition, written without a field"
    assert got[1]["node"] == "pitchfork"
    assert got[2]["mechanism"] == "one-sided write (fold)"
    assert got[4]["mechanism"] == "write to a distant state (subcritical pitchfork)"


def test_a_change_is_drawn_against_the_realization_it_came_from(data_file, tmp_path):
    got = scenario(data_file, [{"do": "look"}, {"do": "edge", "edge": "pf_bias"}, {"do": "undo"},
                               {"do": "param", "name": "eps", "value": 0.3},
                               {"do": "start", "node": "two_spins"}, {"do": "edge", "edge": "spins_field_moved"},
                               {"do": "start", "node": "spin1_easy_axis"},
                               {"do": "edge", "edge": "spin1_classical_aniso"}], tmp_path)
    assert got[0]["reference"] is None                      # nothing was changed yet
    bias = got[1]["reference"]                              # the landscape and the states along the control before the bias
    assert bias["label"].startswith("pitchfork normal form") and bias["potential"] and bias["scan"] and bias["states"] == 2
    assert "before the change" in got[1]["legend"]
    assert got[2]["reference"] is None                      # undo returns to the realization itself
    assert got[3]["reference"]["label"].startswith("the listed values") and got[3]["reference"]["potential"]
    moved = got[5]["reference"]                             # the field moved to the second spin: the algebra is larger,
    assert got[5]["lead"] == "algebra larger than su(2)"    # and the signal is drawn against that of the rotation it left
    assert moved["signal"] and moved["label"].startswith("two coupled spins")
    assert "frequenciesofR1:2.236" in got[5]["facts"].replace(" ", "")   # one frequency: the law of X0 is kept
    assert got[7]["reference"] is None                      # a spin and a classical direction have no common plot


def test_the_two_spins_of_chapter_11_rotate(data_file, tmp_path):
    first = scenario(data_file, [{"do": "start", "node": "two_spins"}], tmp_path)[0]
    assert first["node"] == "two_spins" and first["path"] == ["two_spins"]
    assert first["lead"] == "Bloch rotation" and "2.236" in first["facts"]


def test_an_edge_changes_the_realization_and_undo_returns(data_file, tmp_path):
    got = scenario(data_file, [{"do": "start", "node": "two_spins"}, {"do": "edge", "edge": "spins_field_moved"},
                               {"do": "undo"}], tmp_path)
    assert got[1]["node"] == "two_spins_x1" and got[1]["lead"] == "algebra larger than su(2)"
    assert "dimensionofthealgebra6" in got[1]["facts"].replace(" ", "")
    assert got[2]["node"] == "two_spins" and got[2]["path"] == ["two_spins"]


def test_the_sequence_from_the_spins_to_the_magnet_applies_its_steps_in_order(data_file, tmp_path):
    steps = [{"do": "sequence", "id": "spin-to-magnet"}] + [{"do": "next"}] * 12
    got = scenario(data_file, steps, tmp_path)
    nodes = [g["node"] for g in got]
    assert nodes[0] == "two_spins" and nodes[-1] == "sw_easy"
    assert "spin1_easy_axis" in nodes and "stoner_wohlfarth" in nodes
    at = {g["node"]: g for g in got}
    assert at["spin1_easy_axis"]["lead"] == "algebra larger than su(2)"
    assert at["stoner_wohlfarth"]["lead"] == "2 stable states"
    assert at["two_spins_z0"]["lead"] == "conserved observable"
    assert got[-1]["sequence"].startswith("13 / 13")
    back = scenario(data_file, steps + [{"do": "prev"}], tmp_path)[-1]
    assert back["node"] == "sw_oblique" and back["sequence"].startswith("12 / 13")


def test_the_sequence_from_the_normal_form_reaches_seven_mechanisms(data_file, tmp_path):
    got = scenario(data_file, [{"do": "sequence", "id": "canonical"}] + [{"do": "next"}] * 11, tmp_path)
    assert [g["node"] for g in got] == ["pitchfork", "pitchfork_below", "pitchfork", "pitchfork_bias", "pitchfork",
                                         "pitchfork_subcritical", "pitchfork", "stoner_wohlfarth", "spin1_easy_axis",
                                         "spin1_transverse", "two_spins_h0", "two_spins_z0"]
    assert got[-1]["count"].startswith("11 changes · 7 mechanisms")
    assert got[-1]["mechanism"] == "conserved observable"


def test_a_realization_is_reached_by_the_shortest_path_of_changes(data_file, tmp_path):
    got = scenario(data_file, [{"do": "start", "node": "two_spins"}, {"do": "walk", "node": "sw_easy"}], tmp_path)[-1]
    assert got["path"] == ["two_spins", "two_spins_h0", "spin1_transverse", "spin1_easy_axis", "stoner_wohlfarth",
                           "sw_oblique", "sw_easy"]
    assert got["count"].startswith("6 changes")
