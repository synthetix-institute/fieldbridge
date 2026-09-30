"""The data of the web page: every edit changes only the component it names, every text is filled from a
calculation, and the recorded law constants belong to the current specifications."""
import json
import re
from pathlib import Path

import pytest

pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sympy")

from fieldbridge import site_data as sd  # noqa: E402
from fieldbridge import site_registry as reg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FACTS = {
    "unitary": {"status", "word", "dim", "closure", "frequencies", "frequencies_text", "frequency_list", "carrier", "hilbert", "sector",
                "rate", "theta", "rep", "residual",
                "cause", "couplings", "tunnelling", "energy_difference", "pairing_amplitude", "single_particle_energy",
                "field_x", "field_z", "rabi_frequency", "detuning", "ising_coupling_g", "transverse_field_h"},
    "dissipative": {"states", "states_text", "loss", "write_point", "write_param", "write_kind", "scale_control",
                    "operating_param", "operating_value"}
                   | {f"{p}_{k}" for p in ("sym", "thr", "lock")
                      for k in ("status", "word", "class", "obstruction", "law", "law_err")}
                   | {"lock_ratio", "lock_K", "lock_period", "lock_width", "lock_width_err"},
    "field": {"law", "exponent", "exponent_value", "exponent_fit", "rate"},
    "stochastic": {"convention", "growth", "correction"},
}


def family_of(nid):
    d = next(n for n in reg.NODES if n["id"] == nid)
    if "family" in d:
        return d["family"]
    return "unitary" if "attach" in d else family_of(d["base"])


# ------------------------------------------------------------------------------------------------ the registry
def test_every_edge_and_step_refers_to_the_registry():
    ids = {n["id"] for n in reg.NODES}
    edges = {e["id"]: e for e in reg.EDGES}
    assert len(edges) == len(reg.EDGES)
    assert all(e["from"] in ids and e["to"] in ids and e["slot"] in reg.SLOTS for e in reg.EDGES)
    for s in reg.SEQUENCES:
        for st in s["steps"]:
            assert ("node" in st and st["node"] in ids) or ("edge" in st and st["edge"] in edges) or "prepare" in st
    assert set(reg.SHORT) == ids


def test_every_placeholder_names_a_computed_fact():
    for e in reg.EDGES:
        for m in sd.PLACEHOLDER.finditer(e["text"]):
            nid = e["from"] if m.group(1) else e["to"]
            assert m.group(2) in FACTS[family_of(nid)], (e["id"], m.group(0))
    for s in reg.SEQUENCES:
        at = None
        for st in s["steps"]:
            at = st.get("node") or (next(e for e in reg.EDGES if e["id"] == st["edge"])[
                "from" if st.get("reverse") else "to"] if "edge" in st else at)
            for m in sd.PLACEHOLDER.finditer(st.get("text", "")):
                assert m.group(2) in FACTS[family_of(at)], (s["id"], m.group(0))


def test_every_tutorial_link_points_to_a_section():
    def slug(heading):
        h = heading.strip().lower()
        h = re.sub(r"[^\w\- ]", "", h)
        return h.replace(" ", "-")
    links = [n["tutorial"] for n in reg.NODES if n.get("tutorial")]
    links += [p["tutorial"] for p in reg.PLANNED if p.get("tutorial")] + [reg.M4_RETENTION]
    for link in links:
        path, _, anchor = link.partition("#")
        text = (ROOT / path).read_text(encoding="utf-8")
        if anchor:
            slugs = {slug(line.lstrip("#")) for line in text.splitlines() if line.startswith("#")}
            assert anchor in slugs, (link, anchor)


def test_texts_are_filled_from_facts_and_a_missing_fact_is_an_error():
    assert sd.fill("rate {rate}, angle {theta}°", {"rate": 2.2360679, "theta": 63.4349}) == "rate 2.24, angle 63.4°"
    assert sd.fill("{exponent} and {from.states_text}", {"exponent": "−1/2"}, {"states_text": "two"}) == "−1/2 and two"
    with pytest.raises(sd.SiteError):
        sd.fill("{rate}", {})
    assert sd.fill("{rate}", {}, strict=False) == "—"


# ------------------------------------------------------------------------------------------------ the edits
NON_CODISCOVERY = [e for e in reg.EDGES if e["kind"] != "codiscovery"]


@pytest.fixture(scope="module")
def nodes():
    return sd.resolve()


@pytest.mark.parametrize("edge", NON_CODISCOVERY, ids=[e["id"] for e in NON_CODISCOVERY])
def test_every_edit_changes_only_the_component_it_names(edge, nodes):
    assert sd.check_edit(edge, nodes, {})


def test_an_edit_that_changes_another_component_is_rejected(nodes):
    base = next(e for e in reg.EDGES if e["id"] == "toggle_promoters")
    for wrong in ({**base, "reduces": None}, {**base, "reduces": {"gamma": 1.1}}, {**base, "to": "repressor_ring4"},
                  {**base, "kind": "term", "slot": "Omega", "to": "repressor_ring4"}, {**base, "slot": "Omega"}):
        with pytest.raises(sd.SiteError):
            sd.check_edit(wrong, nodes, {})
    field = next(e for e in reg.EDGES if e["id"] == "field_conservation")
    with pytest.raises(sd.SiteError):
        sd.check_edit({**field, "to": "field_dipole_1d"}, nodes, {})  # conservation and the write at once
    spin = next(e for e in reg.EDGES if e["id"] == "spins_measure_z")
    with pytest.raises(sd.SiteError):
        sd.check_edit({**spin, "to": "two_spins_x1"}, nodes, {})  # the Hamiltonian changes too


# ------------------------------------------------------------------------------------------------ the record
def test_recorded_law_constants_belong_to_the_current_specifications(nodes):
    record = sd.read_law_record()
    assert record, "docs/site/law_constants.json is missing: python3 -B -m fieldbridge demo --law --save-law-record"
    for nid, entry in record["nodes"].items():
        assert nid in nodes, nid
        assert entry["spec_sha256"] == sd.spec_hash(nodes[nid]["spec"]), (
            f"{nid}: the specification changed since docs/site/law_constants.json was written; regenerate it with "
            "python3 -B -m fieldbridge demo --law --save-law-record")
        for target, law in entry.items():
            if target != "spec_sha256":
                assert law["stderr"] >= 0 and abs(law["constant"] - law["expected"]) < 5 * max(law["stderr"], 0.01)


def test_the_canonical_forms_are_the_mechanisms_the_page_names():
    only = ["pitchfork", "pitchfork_below", "pitchfork_bias", "pitchfork_subcritical", "two_spins",
            "rotation_canonical", "rotation_axis"]
    data = sd.build(only=only, derive=False, strict=False, log=lambda *a: None)
    nodes = data["nodes"]
    assert data["start"] == "pitchfork"
    assert {n: nodes[n]["class"] for n in only} == {
        "pitchfork": "symmetric-write", "pitchfork_below": "single-state", "pitchfork_bias": "threshold-write",
        "pitchfork_subcritical": "subcritical-write", "two_spins": "rotation", "rotation_canonical": "rotation",
        "rotation_axis": "conserved"}
    assert all(nodes[n]["universal"] for n in only if n != "two_spins") and not nodes["two_spins"]["universal"]
    # the constant bias h unfolds the pitchfork into a fold at eps = 3 (h/2)^(2/3); the subcritical form jumps at 0
    assert nodes["pitchfork_bias"]["facts"]["write_point"] == pytest.approx(3 * 0.1 ** (2 / 3), abs=1e-3)
    assert nodes["pitchfork_subcritical"]["facts"]["write_kind"] == "subcritical pitchfork"
    assert abs(nodes["pitchfork_subcritical"]["facts"]["write_point"]) < 1e-3
    # the canonical rotation keeps the rate and the angle of the two spins it was detached from
    for k in ("rate", "theta"):
        assert nodes["rotation_canonical"]["facts"][k] == pytest.approx(nodes["two_spins"]["facts"][k], rel=1e-9)
    assert {k: m["node"] for k, m in data["mechanisms"].items()} == {
        "rotation": "rotation_canonical", "conserved": "rotation_axis", "single-state": "pitchfork_below",
        "symmetric-write": "pitchfork", "threshold-write": "pitchfork_bias", "subcritical-write": "pitchfork_subcritical"}


def test_memory_in_model_materials_lists_the_materials_of_the_memory_examples():
    """The section on memory lists every specification file of examples/memory that is not a field, in the order of
    the mechanism that writes it; variants built from a material, and the fields, are not materials."""
    only = ["pitchfork", "pitchfork_bias", "toggle", "colloid_patch", "brusselator", "lotka_volterra",
            "field_nonconserved"]
    data = sd.build(only=only, log=lambda *a: None)
    nodes = data["nodes"]
    assert data["materials"] == ["toggle", "pitchfork", "colloid_patch", "brusselator", "lotka_volterra"]  # writes first
    assert [nodes[i]["class"] for i in data["materials"]] == ["symmetric-write", "symmetric-write", "field-write",
                                                             "oscillation", "neutral-cycles"]
    assert nodes["colloid_patch"]["short"] == "capillary rotors"      # the name of the specification, not "rods"
    assert nodes["toggle"]["facts"]["loss"] == "activation between stored states (Law 3)"
    assert "neutral cycles" in nodes["lotka_volterra"]["facts"]["loss"]  # as on the memory card, not a limit cycle
    # every material of the repository, in the full registry
    files = {p.relative_to(sd.ROOT).as_posix() for p in (sd.ROOT / "examples/memory").rglob("*.json")
             if "fields" not in p.parts}
    named = {d["spec"] for d in sd.all_defs() if "spec" in d}
    assert files <= named and len(files) == 21
    assert reg.SEQUENCES[0]["id"] == "memory-writes"                  # the guided sequences start with memory
    # the outcomes without memory or without a rotation are marked, the retention laws explained, the planned listed
    assert data["absent"] == reg.CLASS_ABSENT and set(reg.CLASS_ABSENT) <= set(reg.CLASSES)
    assert [reg.CLASS_ABSENT[k][0] for k in ("single-state", "neutral-cycles", "conserved", "obstructed")] == [
        "no memory", "no memory", "no rotation", "no rotation"]
    assert sorted(data["retention_laws"]["laws"]) == ["1", "2", "3"]
    assert [q["id"] for q in data["planned"]] == ["frustrated-loops", "retention-rewriting", "hopf-onset", "kuramoto",
                                                   "return-point"]
    assert all(q["law"] and q["exists"] and q["missing"] for q in data["planned"])


def test_every_mechanism_opens_on_a_realization_of_its_class():
    assert set(reg.MECHANISMS) == set(reg.CLASSES)
    ids = {n["id"] for n in reg.NODES}
    assert all(m["node"] in ids for m in reg.MECHANISMS.values()) and reg.START in ids


def test_a_small_site_builds_and_its_sequences_are_complete(tmp_path):
    only = ["two_spins", "two_spins_x1", "two_spins_h0", "two_spins_z0", "field_nonconserved", "field_charge_1d",
            "field_dipole_1d", "log_ito", "log_stratonovich"]
    data = sd.build(only=only, log=lambda *a: None)
    assert set(data["nodes"]) == set(only)
    assert {e["id"] for e in data["edges"]} == {"spins_field_off", "spins_measure_z", "spins_field_moved",
                                                 "field_conservation", "field_dipole_write", "stochastic_convention"}
    assert all("—" not in e["text"] for e in data["edges"])
    assert data["nodes"]["two_spins"]["facts"]["rate"] == pytest.approx(5 ** 0.5)
    assert data["nodes"]["two_spins_z0"]["class"] == "conserved"
    assert data["nodes"]["field_dipole_1d"]["facts"]["exponent"] == "−3/2"
    assert [s["id"] for s in data["sequences"]] == []  # every sequence needs nodes outside this subset
    page = sd.write(data, tmp_path / "site" / "data.js")
    assert page.read_text(encoding="utf-8").startswith("window.FIELDBRIDGE_SITE = ")
    from fieldbridge.web_demo import publish_assets
    html = publish_assets(tmp_path).read_text(encoding="utf-8")
    assert "site/expression.js?v=" in html and "site/data.js?v=" in html
    assert "<select" not in html  # no list of models to choose from
