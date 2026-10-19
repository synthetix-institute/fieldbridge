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


def test_material_choices_start_with_realizations_of_the_current_mechanism(data_file, tmp_path):
    before, opened = scenario(data_file, [{"do": "look"}, {"do": "choose-material"}], tmp_path)
    assert before["chooser"]["hidden"]
    chooser = opened["chooser"]
    assert not chooser["hidden"] and chooser["expanded"] == "true"
    assert chooser["mechanism"] == "symmetric-write"
    assert 'data-material-id="stoner_wohlfarth"' in chooser["options"]
    assert 'data-material-id="pitchfork"' not in chooser["options"]
    assert 'data-material-id="two_spins"' not in chooser["options"]
    assert "Source equations" in chooser["options"] and "Derivation" in chooser["options"]


def test_material_choices_include_spins_and_search_all_mechanisms(data_file, tmp_path):
    opened, all_rows, spins, empty = scenario(data_file, [
        {"do": "choose-material"}, {"do": "material-filter", "key": "mechanism", "value": ""},
        {"do": "material-filter", "key": "search", "value": "spins"},
        {"do": "material-filter", "key": "search", "value": "no such material"}], tmp_path)
    assert 'data-material-id="two_spins"' not in opened["chooser"]["options"]
    assert 'data-material-id="two_spins"' in all_rows["chooser"]["options"]
    assert 'data-material-id="two_spins"' in spins["chooser"]["options"]
    assert 'data-material-id="stoner_wohlfarth"' not in spins["chooser"]["options"]
    assert empty["chooser"]["count"] == "0 model realizations"
    assert "No stored realization matches" in empty["chooser"]["options"]


def test_choosing_a_material_loads_it_directly_and_leaves_the_guide(data_file, tmp_path):
    got = scenario(data_file, [{"do": "sequence", "id": "canonical"}, {"do": "choose-material"},
                               {"do": "material-select", "node": "two_spins"}], tmp_path)[-1]
    assert got["node"] == "two_spins" and got["path"] == ["two_spins"]
    assert got["chooser"]["hidden"] and got["chooser"]["expanded"] == "false"


def test_material_field_filter_and_keyboard_close(data_file, tmp_path):
    got = scenario(data_file, [{"do": "choose-material"},
        {"do": "material-filter", "key": "mechanism", "value": ""},
        {"do": "material-filter", "key": "field", "value": "magnetism"},
        {"do": "material-escape"}, {"do": "choose-material"}, {"do": "material-close"}], tmp_path)
    assert 'data-material-id="stoner_wohlfarth"' in got[2]["chooser"]["options"]
    assert 'data-material-id="two_spins"' not in got[2]["chooser"]["options"]
    assert got[3]["chooser"]["hidden"] and not got[4]["chooser"]["hidden"]
    assert got[4]["chooser"]["mechanism"] == "symmetric-write"
    assert got[5]["chooser"]["hidden"]


def test_the_page_names_the_memory_materials_it_contains(data_file, tmp_path):
    first = scenario(data_file, [{"do": "look"}], tmp_path)[0]
    memory = first["memory"]
    assert memory.startswith("2 materials from 2 fields. Retention: 2 activation between stored states (Law 3)")
    for name in ("pitchfork normal form", "Stoner-Wohlfarth particle", "symmetric write (supercritical pitchfork)",
                 "2 materials", "How is a bit written and retained at a supercritical pitchfork?"):
        assert name in memory, name
    # every mention of a retention law links to its definition, and the three laws are explained beside the materials
    assert first["memoryHtml"].count('class="law-link"') == 3 and "18_memory_writing_and_retention.md#13-retention" in first["memoryHtml"]
    assert "Law 1 relaxation in a curved minimum" in first["laws"] and "Law 3 activation over a barrier" in first["laws"]
    assert "stops decreasing, at ½ ln(1 + W)" in first["laws"] and "Law 3 needs wells separated by a barrier" in first["laws"]
    # the map opens on the memory module and hatches its outcome without memory, one state
    assert first["module"] == "memory" and first["absent"] == 1
    # the mechanisms in preparation
    for name in ("frustrated loops", "retention against rewriting", "onset of oscillation (Hopf)", "What is missing"):
        assert name in first["planned"], name


def test_the_map_shows_one_module_and_follows_the_selection(data_file, tmp_path):
    got = scenario(data_file, [{"do": "look"}, {"do": "module", "key": "all"}, {"do": "walk", "node": "two_spins"},
                               {"do": "module", "key": "memory"}, {"do": "walk", "node": "pitchfork"},
                               {"do": "walk", "node": "two_spins"}], tmp_path)
    # the page opens on the memory module; all mechanisms add the quantum cards and their two hatched outcomes
    # (conserved, several frequencies) to the one without memory
    assert got[0]["module"] == "memory" and got[0]["absent"] == 1
    assert got[1]["module"] == "all" and got[1]["absent"] == 3 and got[1]["cards"] > got[0]["cards"]
    # the whole map stays while the selection moves to a quantum mechanism
    assert got[2]["module"] == "all" and got[2]["node"] == "two_spins"
    # a module chosen by hand is kept while the selection is elsewhere, and while it moves within that module
    assert got[3]["module"] == "memory" and got[3]["node"] == "two_spins"
    assert got[4]["module"] == "memory" and got[4]["absent"] == 1
    # the map follows the selection into another module
    assert got[5]["module"] == "quantum" and got[5]["absent"] == 2


def test_one_change_of_the_normal_form_gives_another_mechanism(data_file, tmp_path):
    got = scenario(data_file, [{"do": "edge", "edge": "pf_below"}, {"do": "undo"}, {"do": "edge", "edge": "pf_bias"},
                               {"do": "undo"}, {"do": "edge", "edge": "pf_subcritical"}], tmp_path)
    assert (got[0]["mechanism"], got[0]["lead"]) == ("single stable state", "one stable state")
    # the selection is shown beside the expression at the top of the page
    assert got[0]["selected"] == "single stable state | pitchfork normal form below the transition, written without a field"
    assert got[1]["node"] == "pitchfork"
    assert got[2]["mechanism"] == "one-sided write (fold)"
    assert got[4]["mechanism"] == "write to a distant state (subcritical pitchfork)"


def test_visible_constructor_controls_restore_the_model(data_file, tmp_path):
    got = scenario(data_file, [{"do": "edge", "edge": "pf_below"},
                               {"do": "toolbar-undo"},
                               {"do": "edge", "edge": "pf_bias"},
                               {"do": "reset"}], tmp_path)
    assert got[0]["lead"] == "one stable state"
    assert got[1]["node"] == "pitchfork" and got[1]["lead"] == "2 stable states"
    assert got[3]["path"] == ["pitchfork"]
    assert "Calculated in the browser" in got[3]["facts"]
    assert "Stored Python calculation" in got[3]["facts"]


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
    assert got[5]["lead"] == "Bloch rotation"               # the rotation is reached through the closure of X0,
    assert moved["signal"] and moved["label"].startswith("two coupled spins")   # and drawn against the two spins
    assert "frequenciesofR1:2.236" in got[5]["facts"].replace(" ", "")   # one frequency: the law of X0 is kept
    assert got[7]["reference"] is None                      # a spin and a classical direction have no common plot


def test_the_two_spins_of_chapter_11_rotate(data_file, tmp_path):
    first = scenario(data_file, [{"do": "start", "node": "two_spins"}], tmp_path)[0]
    assert first["node"] == "two_spins" and first["path"] == ["two_spins"]
    assert first["lead"] == "Bloch rotation" and "2.236" in first["facts"]


def test_an_edge_changes_the_realization_and_undo_returns(data_file, tmp_path):
    got = scenario(data_file, [{"do": "start", "node": "two_spins"}, {"do": "edge", "edge": "spins_field_moved"},
                               {"do": "undo"}], tmp_path)
    assert got[1]["node"] == "two_spins_x1" and got[1]["lead"] == "Bloch rotation"
    assert "dimensionofthealgebra6" in got[1]["facts"].replace(" ", "")
    assert got[2]["node"] == "two_spins" and got[2]["path"] == ["two_spins"]


def test_the_rotation_is_reached_through_the_closure_of_the_observable(data_file, tmp_path):
    """With the field on the second spin, H and X0 generate an algebra of dimension 6 while the closure of X0 keeps
    three operators: the page reports the rotation, its derivation has the letter O, and the axis drawn on the sphere
    is that of the two spins. With the field on both spins the closure has five operators; switching the field on
    the first spin off in the browser gives the first case again."""
    got = scenario(data_file, [{"do": "start", "node": "two_spins"}, {"do": "start", "node": "two_spins_x1"},
                               {"do": "start", "node": "two_spins_both"}, {"do": "term", "k": 1, "on": False}], tmp_path)
    spins, moved, both, off = got
    assert spins["derivation"].startswith("DerivationAKL") and moved["derivation"].startswith("DerivationAOKL")
    assert moved["lead"] == "Bloch rotation" and "close on three operators" in moved["text"]
    assert "dimension 6, not su(2)" in moved["text"] and "the algebra is su(2)" in moved["text"]   # the term is named
    assert moved["axis"] == pytest.approx(spins["axis"], abs=1e-9)
    assert sum(x * x for x in moved["axis"]) == pytest.approx(5.0)                  # |Omega|^2 = 4 g^2 + 4 h^2
    assert both["lead"] == "several frequencies" and both["derivation"].startswith("DerivationAO·")
    assert "5 operators, which move with 2 frequencies" in both["text"]
    assert off["lead"] == "Bloch rotation" and "close on three operators" in off["text"]
    assert off["axis"] == pytest.approx(spins["axis"], abs=1e-9)


def test_following_a_construction_lists_those_through_the_material_and_leaving_returns_to_it(data_file, tmp_path):
    """The button opens the list of constructions, those through the realization on screen first; a construction
    replaces the visitor's material, and leaving it restores the material, the path of changes and the settings."""
    got = scenario(data_file, [
        {"do": "start", "node": "pitchfork"}, {"do": "edge", "edge": "pf_bias"},
        {"do": "param", "name": "eps", "value": 0.3},
        {"do": "open-sequences"},
        {"do": "sequence", "id": "spin-to-magnet"}, {"do": "next"}, {"do": "next"},
        {"do": "leave"}], tmp_path)
    opened = got[3]["sequences"]
    assert opened["expanded"] == "true" and not opened["menuHidden"]
    first = opened["picker"].split("</button>")[0]
    assert 'data-seq="canonical"' in first and "passes through" in first
    during = got[6]
    assert during["node"] != "pitchfork_bias" and during["sequences"]["title"].startswith("From a spin rotation")
    assert during["sequences"]["expanded"] == "false" and not during["sequences"]["playerHidden"]
    back = got[-1]
    assert back["node"] == "pitchfork_bias" and back["path"] == ["pitchfork", "pitchfork_bias"]
    assert back["params"]["eps"] == 0.3 and back["sequences"]["playerHidden"]


def test_the_guide_opens_on_a_first_visit_and_a_realization_from_a_list_opens_directly_with_a_way_back(data_file, tmp_path):
    got = scenario(data_file, [{"do": "look"}, {"do": "guide-close"}, {"do": "guide-toggle"},
                               {"do": "edge", "edge": "pf_bias"}, {"do": "open", "node": "two_spins"}, {"do": "back"}],
                   tmp_path)
    first, closed, reopened, _, opened, back = got
    assert not first["guide"]["hidden"] and first["guide"]["expanded"] == "true"
    assert closed["guide"]["hidden"] and closed["guide"]["expanded"] == "false"
    assert not reopened["guide"]["hidden"]
    # a realization chosen in a list is shown directly, not reached by walking the chain of changes
    assert opened["node"] == "two_spins" and opened["path"] == ["two_spins"]
    assert not opened["back"]["hidden"] and opened["back"]["text"].startswith("↩ Back to")
    assert back["back"]["hidden"]


def test_the_sequence_from_the_spins_to_the_magnet_applies_its_steps_in_order(data_file, tmp_path):
    steps = [{"do": "sequence", "id": "spin-to-magnet"}] + [{"do": "next"}] * 12
    got = scenario(data_file, steps, tmp_path)
    nodes = [g["node"] for g in got]
    assert nodes[0] == "two_spins" and nodes[-1] == "sw_easy"
    assert "spin1_easy_axis" in nodes and "stoner_wohlfarth" in nodes
    at = {g["node"]: g for g in got}
    assert at["spin1_easy_axis"]["lead"] == "several frequencies"
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


def test_hysterons_draw_a_loop_and_a_change_of_couplings_breaks_the_return(tmp_path, monkeypatch):
    monkeypatch.setattr(sd, "HYSTERON_SUBLOOPS", 4)
    data = sd.build(only=["rfim_ferromagnet", "rfim_antiferromagnet", "rfim_antiferromagnet_staggered"], derive=False,
                    strict=False, log=lambda *a: None)
    small = sd.write(data, tmp_path / "site" / "data.js")
    got = scenario(small, [{"do": "start", "node": "rfim_ferromagnet"}, {"do": "edge", "edge": "hyst_couplings_reversed"},
                           {"do": "edge", "edge": "hyst_drive_staggered"}], tmp_path)
    assert [g["lead"] for g in got] == ["return-point memory", "return not exact", "return-point memory"]
    assert "the return is guaranteed" in got[0]["facts"] and "512 of 512 couplings" in got[1]["facts"]
    assert "major loop" in got[0]["legend"] and "unlike the state at H₁" in got[1]["legend"]
    assert got[1]["reference"] and "Grey, dashed" in got[1]["legend"]   # the ferromagnet's loop, before the change
    assert "Relabeling one sublattice" in got[2]["text"]


def test_a_population_that_leaves_its_state_after_a_step(tmp_path):
    data = sd.build(only=["dec_cim_write", "dec_cim_step"], derive=False, strict=False, log=lambda *a: None)
    data_file = sd.write(data, tmp_path / "data.js")
    got = scenario(data_file, [{"do": "start", "node": "dec_cim_write"}, {"do": "walk", "node": "dec_cim_step"}],
                   tmp_path)
    assert got[0]["module"] == "decision" and got[0]["mechanism"] == "follows the bias (swept write)"
    assert "along the actual passage" in got[0]["legend"] and "threshold of the symmetric state" in got[0]["facts"]
    assert got[1]["node"] == "dec_cim_step" and got[1]["mechanism"] == "passage from a one-component seed"
    assert "ln 13.09 = 2.5718" in got[1]["facts"] and "pump ramp → pump step" in got[1]["changes"] + got[0]["changes"]
