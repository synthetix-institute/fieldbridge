from __future__ import annotations

from .models import AnalogyMatch, Comparison, Fingerprint, MechanismRecord, MechanismSheet, Translation
from .routes import FIBERS, ROUTES


def is_hyperion_record(record: MechanismRecord) -> bool:
    return record.field_id == "hyperion_equation" or record.source == "hyperion_equation_witnesses"


def hyperion_reference(record: MechanismRecord) -> str:
    hyperion = record.hyperion or {}
    paper_id = hyperion.get("paper_id") or next((ref.replace("arxiv:", "") for ref in record.references if ref.startswith("arxiv:")), "")
    arxiv_url = hyperion.get("arxiv_url") or (f"https://arxiv.org/abs/{paper_id}" if paper_id else "")
    apparatus = hyperion.get("apparatus_regime") or "Hyperion witness"
    omega = hyperion.get("omega_tokens") or ""
    witness_id = hyperion.get("witness_id") or record.record_id
    link = f"[{paper_id}]({arxiv_url})" if paper_id and arxiv_url else "source unavailable"
    activation = f"{apparatus} {omega}".strip()
    return f"{witness_id}: {activation}; arXiv {link}"


def render_equation_rows(record: MechanismRecord) -> list[str]:
    if not record.equations:
        if is_hyperion_record(record):
            return ["- raw equation snippet hidden; use the arXiv audit link for source context"]
        return ["- no equations supplied"]
    return [f"- `{equation}`" for equation in record.equations]


def render_fingerprint(fp: Fingerprint) -> str:
    rows = ["# Fingerprint", "", "## Routes"]
    for route_id, score in sorted(fp.routes.items(), key=lambda item: item[1], reverse=True):
        rows.append(f"- `{route_id}` {score:.2f}: {ROUTES[route_id][0]}")
    rows.extend(["", "## Fibers"])
    for fiber_id, score in sorted(fp.fibers.items(), key=lambda item: item[1], reverse=True):
        rows.append(f"- `{fiber_id}` {score:.2f}: {FIBERS[fiber_id][0]}")
    return "\n".join(rows)


def render_match(match: AnalogyMatch) -> str:
    record = match.record
    routes = ", ".join(f"{key}:{value:.2f}" for key, value in match.route_overlap.items()) or "none"
    fibers = ", ".join(f"{key}:{value:.2f}" for key, value in match.fiber_overlap.items()) or "none"
    hits = ", ".join(match.keyword_hits) or "none"
    return "\n".join(
        [
            f"## {record.title}",
            "",
            f"- score: `{match.score:.3f}`",
            f"- field: `{record.field_id}`",
            f"- route overlap: {routes}",
            f"- fiber overlap: {fibers}",
            f"- keyword hits: {hits}",
            f"- summary: {record.summary}",
            *( [f"- audit link: {hyperion_reference(record)}"] if is_hyperion_record(record) else [] ),
            "",
            "Equations:",
            *render_equation_rows(record),
            "",
            "Controls:",
            *(f"- {control}" for control in record.controls),
        ]
    )


def render_search(matches: list[AnalogyMatch]) -> str:
    if not matches:
        return "No analogs found."
    return "\n\n".join(render_match(match) for match in matches)


def render_mechanism_sheet(sheet: MechanismSheet) -> str:
    rows = [
        "# FieldBridge Mechanism Sheet",
        "",
        f"Source: **{sheet.source_title}**",
        "",
        "## Extracted Invariant",
        "",
        sheet.invariant,
        "",
        "## Roles",
        "",
        f"- state: {sheet.state}",
        f"- input: {sheet.input_signal}",
        f"- boundary: {sheet.boundary}",
        f"- output: {sheet.output}",
        "",
        "## Equations",
        "",
        *(f"- `{item}`" for item in sheet.equations),
        "",
        "## Measurements",
        "",
        *(f"- {item}" for item in sheet.measurements),
        "",
        "## Falsifying Controls",
        "",
        *(f"- {item}" for item in sheet.controls),
        "",
        "## Active Routes",
        "",
        *(f"- {item}" for item in sheet.route_explanation),
        "",
        "## Evidence Boundary",
        "",
        sheet.evidence_boundary,
    ]
    return "\n".join(rows)


def render_comparison(comparison: Comparison) -> str:
    def render_scores(rows: dict[str, float]) -> str:
        return ", ".join(f"{key}:{value:.2f}" for key, value in rows.items()) or "none"

    return "\n".join(
        [
            "# FieldBridge Paper-to-Paper Comparison",
            "",
            f"Similarity score: `{comparison.score:.3f}`",
            "",
            comparison.interpretation,
            "",
            "## Preserved Structure",
            "",
            f"- routes: {render_scores(comparison.preserved_routes)}",
            f"- fibers: {render_scores(comparison.preserved_fibers)}",
            "",
            "## Changed Structure",
            "",
            f"- routes: {render_scores(comparison.changed_routes)}",
            f"- fibers: {render_scores(comparison.changed_fibers)}",
            "",
            "## Source Mechanism",
            "",
            f"- state: {comparison.source.state}",
            f"- boundary: {comparison.source.boundary}",
            f"- output: {comparison.source.output}",
            "",
            "## Target Mechanism",
            "",
            f"- state: {comparison.target.state}",
            f"- boundary: {comparison.target.boundary}",
            f"- output: {comparison.target.output}",
            "",
            "## Analogous Equations",
            "",
            *(f"- `{item}`" for item in comparison.analogous_equations or ["No equation-level analogue crossed the threshold."]),
        ]
    )


def render_translation(translation: Translation) -> str:
    rows = [
        "# FieldBridge Mechanism Translation",
        "",
        f"Target field: **{translation.target_field.label}**",
        "",
        "## Invariant",
        "",
        translation.invariant,
        "",
        "## Target Formulation",
        "",
        translation.target_formulation,
        "",
        "## Variables",
        "",
        *(f"- {item}" for item in translation.variables),
        "",
        "## Mechanism Equations",
        "",
        *(f"- `{item}`" for item in translation.equations),
        "",
        "## Measurements",
        "",
        *(f"- {item}" for item in translation.measurements),
        "",
        "## Falsifying Controls",
        "",
        *(f"- {item}" for item in translation.controls),
        "",
        "## Closest Analogues",
        "",
    ]
    if translation.matches:
        for match in translation.matches:
            if is_hyperion_record(match.record):
                rows.append(f"- `{match.score:.3f}` atlas audit: {hyperion_reference(match.record)}")
            else:
                rows.append(f"- `{match.score:.3f}` field anchor: {match.record.title} ({match.record.field_id})")
    else:
        rows.append("- No field-specific analogue found; rendered from field pack defaults.")
    rows.extend(["", "## Evidence Boundary", "", translation.evidence_boundary])
    return "\n".join(rows)
