# Source Lock

This public repository package is self-contained for replay. All paths below are relative to the repository root.

## Public Sources

| Source | Path | Role |
|---|---|---|
| Manuscript | `paper/main.tex` | Main theorem, type-D parity sublayer, S12 separator, finite appendix |
| Reproduce commands | `REPRODUCE.md` | Replay entrypoint |
| Core scripts | `scripts/p15_s*.py` | Exact theorem certificates |
| Type-D scripts | `type_d/` | Type-D parity sublayer replay scripts + JSON certificates (`opus_g15/g16/g20/g22`); `opus_g22` imports the engine from `scripts/` |
| Appendix scripts | `scripts/p15_full_bn_*.py` | Optional finite `n=8` appendix replay |
| Certified artifacts | `artifacts/certified/` | JSON/Markdown certificate outputs |
| External n=8 atlas inputs | `experiments/external_agents_2026_07_07/` and `artifacts/external_n8_atlas/` | Finite appendix input CSV/metadata |

Private PAPP master maps and exploratory transcripts are not public source dependencies for this repository package.
