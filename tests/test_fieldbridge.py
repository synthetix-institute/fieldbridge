from fieldbridge.routes import fingerprint_text
from fieldbridge.database import load_field_evidence, load_static_index_records
from fieldbridge.extract import compare_mechanisms, extract_mechanism
from fieldbridge.search import find_analogs, translate_mechanism


def test_fingerprint_detects_transport_and_boundary():
    fp = fingerprint_text("partial_t q + nabla dot J = S with boundary B and closure C(q)=0")
    assert fp.routes["transport_flow_route"] > 0.2
    assert fp.routes["boundary_weak_form_route"] > 0.2
    assert fp.routes["constraint_closure_route"] > 0.2


def test_search_returns_material_analog():
    text = "A retained tissue state q changes later repair after washout through boundary closure."
    matches = find_analogs(text, target_field="material_intelligence")
    assert matches
    assert matches[0].record.field_id == "material_intelligence"


def test_translate_has_controls():
    text = "A prior signal writes q; after erasure, y changes through C(q,u,B)=0."
    translation = translate_mechanism(text, target_field="collective_intelligence")
    assert translation.controls
    assert translation.equations
    assert translation.target_field.field_id == "collective_intelligence"


def test_hyperion_static_index_can_support_search():
    records = load_static_index_records(max_records=3)
    assert records
    assert records[0].field_id == "hyperion_equation"
    matches = find_analogs(
        "boundary closure transport equation",
        target_field="material_intelligence",
        include_hyperion=True,
    )
    assert matches


def test_field_evidence_sidecars_load():
    evidence = load_field_evidence()
    assert "material_intelligence" in evidence
    assert evidence["material_intelligence"]["telemetry"]["atlas_candidate_witnesses"] > 0


def test_extract_returns_mechanism_sheet():
    text = "partial_t q = W(u,t)-q/tau; C(q,u,B)=0; y=A(q,B). A wound boundary changes later repair."
    sheet = extract_mechanism(text)
    assert sheet.equations
    assert "state" in sheet.state.lower() or "q" in sheet.state.lower()
    assert sheet.controls


def test_compare_preserves_routes():
    source = "partial_t q = W(u,t)-q/tau; C(q,u,B)=0; y=A(q,B)."
    target = "partial_t m = W(T,t)-m/tau; C(m,u,B)=0; y=A(m,B)."
    comparison = compare_mechanisms(source, target)
    assert comparison.score > 0.4
    assert comparison.preserved_routes
