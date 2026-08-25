from fieldbridge import cli
from fieldbridge.constructor import construct_transfer
from fieldbridge.continuation import validate_future_state
from fieldbridge.routes import fingerprint_text
from fieldbridge.database import load_field_evidence, load_static_index_records
from fieldbridge.extract import compare_mechanisms, extract_mechanism
from fieldbridge.search import find_analogs, translate_mechanism
from fieldbridge.zero_shot import validate_full_paper_zero_shot


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
    assert translation.protocols
    assert translation.equations
    assert translation.target_field.field_id == "collective_intelligence"


def test_cli_pdf_input_uses_shared_document_parser(tmp_path, monkeypatch):
    path = tmp_path / "paper.pdf"
    path.write_bytes(b"%PDF-placeholder")
    monkeypatch.setattr(cli, "read_document", lambda value: f"parsed:{value.name}")

    assert cli.read_input(str(path)) == "parsed:paper.pdf"


def test_constructor_exposes_identity_completion_and_falsifier():
    text = "partial_t p = div(p grad U) + T Laplacian p; integral p = 1; p_inf = exp(-U/T)."
    transfer = construct_transfer(text, target_field="stochastic_optimization", include_hyperion=False)
    assert transfer.readiness == "structurally_complete_constructor_proposal"
    assert transfer.source_identity["M"]["Omega_proxy"]
    assert transfer.target_identity["M"]["Xi_proxy"]
    assert transfer.required_attachments["C_closure"]
    assert transfer.required_attachments["R_readout"]
    assert transfer.required_attachments["P_protocol"]
    assert transfer.required_attachments["falsifiers"]
    assert transfer.target_identity["F"]["P"] == transfer.required_attachments["P_protocol"]
    assert set(transfer.required_attachments["P_protocol"]).isdisjoint(
        transfer.required_attachments["falsifiers"]
    )
    assert {move["acts_on"] for move in transfer.constructor_moves} >= {"Omega", "Xi", "C", "R", "P"}
    assert transfer.validation_gates["independent_target_validation"] is False


def test_tex_display_math_is_not_a_commutator():
    fp = fingerprint_text(r"\[p_\infty \propto \exp[-U/T]\] with probability conservation")
    assert fp.routes["commutator_incompatibility_route"] == 0.0


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


def test_full_paper_zero_shot_excludes_query_paper(tmp_path):
    papers = [
        ("continuum-a", "continuum", "diffusion", "partial_t u = kappa nabla^2 u. No flux boundary conserves mass."),
        ("graph-a", "graphs", "diffusion", "dot u = -kappa L_G u. Graph Laplacian closure conserves total state."),
        ("continuum-b", "continuum", "oscillation", "partial_t^2 q + omega^2 q = 0. Periodic oscillation."),
        ("network-b", "networks", "oscillation", "dot z = A z. Complex eigenmodes produce network oscillation."),
    ]
    manifest = []
    for paper_id, field_id, mechanism_id, text in papers:
        path = tmp_path / f"{paper_id}.txt"
        path.write_text((text + "\n\n") * 30, encoding="utf-8")
        manifest.append({
            "paper_id": paper_id,
            "path": path.name,
            "field_id": field_id,
            "mechanism_id": mechanism_id,
        })
    manifest_path = tmp_path / "manifest.json"
    import json
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    report = validate_full_paper_zero_shot(manifest_path, top_k=2, bootstrap_samples=50)
    assert report["readiness"] == "usable"
    assert report["eligible_queries"] == 4
    assert report["protocol"]["lexical_fit"].startswith("gallery_only")
    assert report["performance_claim_gate"]["passed"] is False
    for query in report["queries"]:
        assert all(item["paper_id"] != query["paper_id"] for item in query["operational_ranking"])
        assert any(item["relevant"] for item in query["operational_ranking"])


def test_future_state_validation_predicts_withheld_equations_not_current_labels(tmp_path):
    rows = []
    for paper_index in range(80):
        first_move = f"T{paper_index % 2}"
        for step in range(5):
            rows.append({
                "paper_id": f"paper-{paper_index:03d}",
                "current_omega": f"Omega{paper_index % 3}",
                "current_xi": f"Xi{step % 2}",
                "current_completion": "C+R" if step % 2 else "C",
                "first_move": first_move,
                "next_move": f"T{(paper_index + step) % 3}",
                "destination_omega": f"Omega{paper_index % 3}",
                "destination_xi": f"Xi{(step + 1) % 2}",
                "destination_completion": "C+R+P" if step % 2 else "C+R",
            })
    manifest = tmp_path / "transitions.jsonl"
    import json
    manifest.write_text(
        "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8"
    )

    report = validate_future_state(manifest, min_evaluation_transitions=5)
    assert report["readiness"] == "usable"
    assert report["protocol"]["unit_of_holdout"] == "complete_paper"
    assert report["protocol"]["current_state_classification"] is False
    assert report["evaluation_papers"] > 0
    assert report["targets"]["destination_omega"]["accuracy"] > 0.9


def _example_text():
    from pathlib import Path

    return Path("examples/brownian_probability_flow.tex").read_text()


def test_atlas_witnesses_reachable_from_every_field_pack():
    """Atlas witnesses are cross-field evidence and must not be hidden by target_fields.

    The shipped index enumerates only the three intelligence fields, so gating
    witnesses on target_fields made the whole 2,633-record index invisible to
    stochastic_optimization -- the field the README's own example targets.
    """
    from fieldbridge.database import load_all
    from fieldbridge.search import is_hyperion_record

    packs, _ = load_all(None)
    for pack in packs:
        matches = find_analogs(
            _example_text(), target_field=pack.field_id, top_k=8, include_hyperion=True
        )
        assert any(is_hyperion_record(m.record) for m in matches), pack.field_id


def test_include_hyperion_changes_the_constructor_result():
    """The flag must alter the transfer; it silently did nothing before."""
    with_atlas = construct_transfer(_example_text(), "stochastic_optimization")
    without = construct_transfer(
        _example_text(), "stochastic_optimization", include_hyperion=False
    )
    assert with_atlas.atlas is not None
    assert without.atlas is None
    assert with_atlas.source_identity["M"]["Omega"] != without.source_identity["M"]["Omega"]


def test_constructor_reports_atlas_tokens_not_only_proxies():
    transfer = construct_transfer(_example_text(), "stochastic_optimization")
    assert "Ω" in transfer.source_identity["M"]["Omega"]
    assert transfer.atlas["witness_id"].startswith("EW")


def test_target_carrier_token_is_not_inherited_from_the_source_witness():
    """Reattachment proposes a target carrier; the atlas never assigned one."""
    transfer = construct_transfer(_example_text(), "stochastic_optimization")
    assert transfer.source_identity["M"]["Xi_token"]
    assert transfer.target_identity["M"]["Xi_token"] is None


def test_constructor_renders_retrieved_witnesses():
    from fieldbridge.render import render_constructor

    out = render_constructor(construct_transfer(_example_text(), "stochastic_optimization"))
    assert "## Retrieved Witnesses" in out
    assert "arXiv" in out
