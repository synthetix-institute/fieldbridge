"""The heredity card of a body that grows and divides: the predictions from its equations beside, on request, lineages
of the full body with noise.

The card holds the threshold and its reduction (L_c, a, b), the conditions of the law (the dip, ln G over a generation,
Lambda), the class, and for a body that inherits through the threshold the law for the daughter of its noise-free
lineage (P in the body, history and normal-form versions). With `simulate` it adds lineages of the full body: the
fraction of daughters that keep their parent's sign against the law averaged over their states at birth, with the gate
|P - P_law| < 4 stderr + 0.02.

Classes (names after the conditions of the law): inherited through a threshold; kept above the threshold (no dip:
the barrier keeps the order); lost in the dip (ln G < 0); threshold moved by division.
"""
from __future__ import annotations

import time
from typing import Dict

from .lineage import agrees, condition
from .predict import predict, reduction


def card(lin, simulate: bool = False, lineages: int = 1000, generations: int = 3, seed: int = 1,
         noise: float = 1.0) -> Dict:
    red = reduction(lin)
    pr = predict(lin, red)
    out = {"name": lin.name, "field": lin.field, "source": lin.source, "path": lin.path,
           "spec_sha256": lin.spec_sha256, "body": lin.kind, "size": lin.size,
           "growth": {"law": lin.growth.law, "rate": lin.growth.rate}, "noise": lin.noise_text,
           "class": pr["class"], "prediction": pr}
    if lin.expect is not None:
        out["expect"] = lin.expect
    if simulate:
        t0 = time.time()
        row = condition(lin, red, pr["L_div"], noise=noise, n=lineages, generations=generations, burn=1, seed=seed)
        row["agrees_body"] = agrees(row, "law_body")
        row["agrees_history"] = agrees(row, "law_history")
        row["seconds"] = time.time() - t0
        out["lineages"] = row
    return out
