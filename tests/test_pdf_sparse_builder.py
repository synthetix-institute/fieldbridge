import json

from fieldbridge.pdf_sparse_builder import build_pdf_field_pack


def test_pdf_sparse_builder_exports_typed_mechanism_kg(tmp_path):
    papers = tmp_path / "papers"
    papers.mkdir()
    (papers / "memory_boundary.txt").write_text(
        """
        A thermal stimulus writes a retained material state q in an interface.
        After washout, the boundary condition B changes the transport response y.
        The mechanism is tested by reset controls, shuffled history controls,
        and a residual C(q,u,B)=0 with partial_t q = W(u,t)-q/tau.
        The preparation protocol applies two pulses before the readout.

        A second optical input writes a polarization memory. The surface boundary
        gates conductance and motion; erasure should remove the later readout.
        """,
        encoding="utf-8",
    )

    summary = build_pdf_field_pack(
        pdf_dir=papers,
        field_id="test_field",
        label="Test Field",
        out_dir=tmp_path / "out",
        max_docs=5,
        max_anchors=4,
        extensions=(".txt",),
    )

    assert summary["documents"] == 1
    evidence_path = tmp_path / "out" / "field_pack_evidence" / "test_field.json"
    graph_path = tmp_path / "out" / "kg" / "test_field_knowledge_graph.json"
    pack_path = tmp_path / "out" / "field_packs" / "test_field.json"
    adapter_path = tmp_path / "out" / "field_adapters" / "test_field.json"
    report_path = tmp_path / "out" / "reports" / "test_field_adapter.md"
    assert evidence_path.exists()
    assert graph_path.exists()
    assert pack_path.exists()
    assert adapter_path.exists()
    assert report_path.exists()

    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    graph = json.loads(graph_path.read_text(encoding="utf-8"))
    adapter = json.loads(adapter_path.read_text(encoding="utf-8"))
    assert evidence["receptor_graph_path"] == "kg/test_field_knowledge_graph.json"
    assert evidence["field_adapter_path"] == "field_adapters/test_field.json"
    assert graph["nodes"]
    assert graph["edges"]
    relationships = {edge["relationship"] for edge in graph["edges"]}
    assert "writes_state" in relationships
    assert "drives_readout" in relationships
    assert graph["statistics"]["nodes"] > 0
    assert "gap_report" in graph
    assert adapter["artifact_type"] == "fieldbridge_field_adapter"
    assert adapter["constructor_roles"]
    assert "state_or_carrier" in adapter["field_native_receivers"]
    assert "protocol_execution" in adapter["field_native_receivers"]
    assert adapter["field_native_receivers"]["protocol_execution"]
    assert "coordinate_domain" in adapter["substrate_profile"]
    assert "route_to_field_adapter" in adapter
    records = json.loads(
        (tmp_path / "out" / "index" / "core_examples.json").read_text(
            encoding="utf-8"
        )
    )
    assert records[0]["references"][0].startswith("txt:")
    assert all(
        "cellular, tissue" not in value
        for record in records
        for value in record["variables"]
    )


def test_pdf_sparse_builder_skips_unreadable_document(tmp_path):
    papers = tmp_path / "papers"
    papers.mkdir()
    (papers / "usable.txt").write_text(
        "A preparation protocol updates state q and measures response y.",
        encoding="utf-8",
    )
    (papers / "broken.pdf").write_bytes(b"not a pdf")

    summary = build_pdf_field_pack(
        pdf_dir=papers,
        field_id="mixed_input",
        label="Mixed Input",
        out_dir=tmp_path / "out",
        max_docs=5,
        max_anchors=2,
        extensions=(".pdf", ".txt"),
    )

    assert summary["documents_discovered"] == 2
    assert summary["documents"] == 1
    assert summary["documents_failed"] == 1
    evidence = json.loads(
        (tmp_path / "out" / "field_pack_evidence" / "mixed_input.json").read_text(
            encoding="utf-8"
        )
    )
    failures = evidence["source_artifacts"]["extraction_failures"]
    assert failures and failures[0]["path"].endswith("broken.pdf")
