"""Keep the inverse-construction formulas legible in GitHub Markdown."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def test_inverse_construction_uses_a_safe_pair_sum():
    source = (ROOT / "docs/tutorial/14_inverse_construction.md").read_text()
    blocks = re.findall(r"```math\n(.*?)\n```", source, re.S)
    assert len(blocks) == 2
    assert r"\sum_{i=0}^{N-2}\sum_{j=i+1}^{N-1}" in blocks[0]
    assert "<" not in blocks[0]
    assert r"\mathrm{Re}\," in blocks[1]
    for block in blocks:
        depth = 0
        for brace in re.findall(r"(?<!\\)[{}]", block):
            depth += 1 if brace == "{" else -1
            assert depth >= 0
        assert depth == 0


def test_constructor_precedes_the_expandable_map_and_catalogues():
    html = (ROOT / "fieldbridge/web/index.html").read_text()
    assert html.index('id="instrument"') < html.index('id="expression"')
    assert 'id="sequence-start"' in html
    assert '<ol class="howto">' not in html
    for section in ("expression", "memory", "planned", "codiscovery", "realizations", "contribute"):
        assert f'<details id="{section}"' in html
    assert "Tutorial and calculations" in html


def test_choose_material_opens_an_inline_list_not_the_memory_chapter():
    html = (ROOT / "fieldbridge/web/index.html").read_text()
    assert 'id="material-choose"' in html
    assert 'aria-controls="material-chooser"' in html
    assert '<a href="#memory">Choose a material</a>' not in html
    assert html.index('id="material-options"') < html.index('class="instrument"')
    for control in ("material-mechanism", "material-field", "material-search"):
        assert f'id="{control}"' in html
