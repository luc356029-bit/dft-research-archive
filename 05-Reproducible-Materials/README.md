# 05 — Reproducible Materials

This directory is the evidence layer for the DFT reproduction project. It follows the actual research sequence: A → B → C1 → supplementary validation.

## Evidence structure
Each stage separates references, inputs, methods and parameters, code, execution records, outputs, and validation. The goal is to make each result traceable to the exact input and procedure that produced it.

## Evidence principles
- The original Google Drive source folders remain unchanged.
- A, B, and C1 describe the main research sequence; Supplementary-Validation contains cross-stage checks, independent comparisons, attribution, and extended analyses.
- Historical plans, failed runs, and superseded variants are not silently promoted to current production evidence.
- A successful script run, a numerically reproducible result, and a scientifically validated conclusion are different claims.

## Public-distribution exclusions
- VASP POTCAR files and embedded copies are not redistributed.
- Third-party publisher PDFs/HTML and non-commercial foundation-model weights are not redistributed.
- The Materials Project API key formerly hard-coded in the source script is not included; the public screening script reads MP_API_KEY from the environment. Rotate any previously exposed key.

## Reproducibility status
This is a curated public evidence set, not a claim that every historical calculation can be rerun from a clean machine with one command. VASP and licensed PAW datasets remain external prerequisites.


## GitHub publication scope

See [PUBLICATION_STATUS.md](PUBLICATION_STATUS.md) for the current published file count and the large raw-output/evidence bundles that remain outside the GitHub tree.
