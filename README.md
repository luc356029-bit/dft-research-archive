# DFT Research Archive

Public, curated record of a computational reproduction project on lithium-alloy screening and Li-ion transport.

## Research sequence

1. **A — High-throughput screening:** retrieve and filter candidate Li alloys; 77 CIF structures are included.
2. **B — Structure and diffusion pathways:** symmetry, inequivalent Li sites, distance shells, periodic connectivity and candidate hops.
3. **C1 — NEB and transport calculations:** structural relaxation, vacancy construction, path setup, NEB evidence and transport analysis.
4. **Supplementary validation:** endpoint checks, independent comparisons, error attribution and extended analyses.

The detailed, ordered materials are under [05-Reproducible-Materials](05-Reproducible-Materials/README.md).

## Evidence and reproducibility

This is a curated public evidence set, not a claim that every historical calculation can be rerun from a clean machine with one command. Each stage distinguishes inputs, methods, code, execution records, outputs and validation.

- Original Google Drive source folders were left unchanged.
- VASP POTCAR files and embedded copies are excluded.
- Third-party publisher PDFs/HTML and non-commercial model weights are not redistributed.
- Materials Project API credentials are not included. The public screening script reads MP_API_KEY from the environment.
- SHA256 manifest and a verification script are included in the project directory.

See [LICENSE_STATUS.md](05-Reproducible-Materials/LICENSE_STATUS.md), [SOURCE_MAP.md](05-Reproducible-Materials/SOURCE_MAP.md), and [POTCAR_NOT_INCLUDED.md](05-Reproducible-Materials/POTCAR_NOT_INCLUDED.md) before reuse.