# RH concept map

Typed research documentation only. **No claim for or against RH.**

This layer maps mathematical objects, criteria, theorems, TFPT structures,
barriers, and open questions. Catalog attempt records remain external and are
referenced by their stable `path`; they are never duplicated as concept nodes.

## Files

- `schema.json`: draft-07 schema and closed relation/status enums.
- `make_seeds.py`: source for the hand-curated seed ontology.
- `seed_nodes.json`, `seed_edges.json`: generated curated seeds.
- `lexicon.json`: generated case-insensitive alias lexicon.
- `extract_links.py`: catalog mention, kill, and co-occurrence extraction.
- `extracted_edges.json`: generated attempt links.
- `build_map.py`: merge, validation, liveness, and viewer export.
- `rh_concept_map.json`: validated map.
- `rhmap.py`: stdlib query/export CLI.
- `gaps_report.json`: generated G1–G6 report.

## Provenance

Every node definition has a `file:line`, DOI, arXiv id, or classical reference.
Every edge has its own source. A `THEOREM` edge requires a classical source or
a named Lean theorem. Corpus-internal results are `CONDITIONAL` or `HEURISTIC`
unless certified; finite certificates never become RH claims. `CO_OCCURS` is
statistical only. Lean equivalence edges certify logical typing, not positivity.

## Add a concept or relation

1. Edit `make_seeds.py`; use a stable kebab-case id and a ≤60-word definition.
2. Give the node and every relation a precise source.
3. Add unambiguous aliases; avoid short tokens such as `RH`, `E8`, or `PRIME`.
4. Regenerate seeds, extraction, map, and gaps.
5. Never hand-edit generated JSON.

## Build

```bash
python3 rh/catalog/map/make_seeds.py
python3 rh/catalog/map/extract_links.py
python3 rh/catalog/map/extract_links.py --llm  # optional, key + budget required
python3 rh/catalog/map/build_map.py
python3 rh/catalog/map/rhmap.py gaps --json
```

Extraction includes `fragments/part_8.json` when present. `--llm` examines at
most 200 `EQUIVALENCE`/`LOAD_BEARING_OPEN` records; its links are explicitly
`source=llm`, `strength=HEURISTIC`, and extraction remains fully offline without
`OPENAI_API_KEY`.

`proposed_additions_index.json` is the review surface for the classical torus
covering-index theorem, the still-open physical index-modular state, and the
Connes cutoff dichotomy. `proposed_additions_bridge2.json` is the review
surface for the Suzuki screw-function / W2–W4 window chain and the Conrey–Li
de Branges note. Both are proposed documentation only and are not a claim
for or against RH.

## Query

```bash
python3 rh/catalog/map/rhmap.py show weil-positivity
python3 rh/catalog/map/rhmap.py neighbors weil-positivity --rel EQUIVALENT_TO
python3 rh/catalog/map/rhmap.py path tfpt-e8-lattice weil-positivity --avoid-killed
python3 rh/catalog/map/rhmap.py equivalents
python3 rh/catalog/map/rhmap.py gaps
python3 rh/catalog/map/rhmap.py stats
python3 rh/catalog/map/rhmap.py export --format graphml /tmp/rh.graphml
```

Add `--json` anywhere for single-line machine output. Human output is capped at
60 lines. Paths are undirected traversals of sourced mathematical relations;
attempt links, co-occurrence, and ontology classification edges are excluded.
