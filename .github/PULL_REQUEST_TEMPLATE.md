## Summary

<!-- One material, mechanism, analysis or correction, in one or two sentences. -->

## Source

<!-- For a material: the paper or book its equations come from. -->

## Checklist

- [ ] `python3 -B -m fieldbridge memory check <file>` passes its required checks (for a material)
- [ ] `docs/materials.md` and `docs/materials.json` are regenerated with `python3 -B -m fieldbridge memory catalog` (for a material)
- [ ] `python3 -B -m pytest -q -p no:cacheprovider` passes
- [ ] Parameter values chosen for the example, rather than taken from the source, are listed under `assumptions`
- [ ] New behaviour comes with a test
