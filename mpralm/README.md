# mpralm v1.16.0

This directory contains R scripts that support MPRA analysis using the `mpralm` package.

## Contents

- `mpralm_Allele.R` runs allele-specific MPRA modeling to test for differential activity between alleles.
- `mpralm_CRE_test_vs_negctrl.R` compares candidate regulatory elements against negative controls to identify active CREs.

## Usage

Both scripts are written as R scripts intended to be run with `Rscript`. Adjust the inputs inside each script to point to the appropriate count matrix and metadata files before executing, e.g.

```bash
Rscript mpralm_Allele.R
Rscript mpralm_CRE_test_vs_negctrl.R
```
Review the top of each script for required package dependencies and expected input formats.
