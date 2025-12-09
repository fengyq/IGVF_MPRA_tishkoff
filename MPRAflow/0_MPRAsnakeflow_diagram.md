# Regulatory SNP Identification Workflow - MPRAsnakeflow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           INPUT DATA PREPARATION                            │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│   MPRA Library      │    │   Sequencing Data   │    │  Variant Definition │
│                     │    │                     │    │                     │
│ • REF sequences     │    │ • DNA-seq (input)   │    │ • ID  REF  ALT      │
│ • ALT sequences     │    │ • RNA-seq (output)  │    │ • chr1:123 G C      │
│ • Barcodes          │    │ • Paired/Single-end │    │ • chr2:456 A T      │
└─────────────────────┘    └─────────────────────┘    └─────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            CRITICAL DECISION POINT                          │
│                           MAPPING STRATEGY SELECTION                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                          ┌───────────┴───────────┐
                          │                       │
                          ▼                       ▼
            ┌─────────────────────┐    ┌─────────────────────┐
            │    EXACT MODE       │    │   ALIGNMENT MODE    │
            │  (RECOMMENDED)      │    │   (BWA/BBMap)       │
            │                     │    │                     │
            │ ✅ 100% specificity │    │ ⚠️  Risk of cross-   │
            │ ✅ Zero contamination│    │    contamination    │
            │ ❌ Lower yield       │    │ ✅ Higher yield      │
            └─────────────────────┘    └─────────────────────┘
                          │                       │
                          └───────────┬───────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            ASSIGNMENT WORKFLOW                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│  Barcode-Oligo      │    │   Read Assignment   │    │   Quality Control   │
│   Mapping           │    │                     │    │                     │
│                     │    │ DNA reads → REF/ALT │    │ • Mapping stats     │
│ BC1 → G-allele      │    │ RNA reads → REF/ALT │    │ • Assignment rates  │
│ BC2 → C-allele      │    │                     │    │ • Cross-contamination│
│ BC3 → G-allele      │    │ EXACT MODE RESULT:  │    │   check             │
└─────────────────────┘    │ G-reads → G-allele  │    └─────────────────────┘
                           │ C-reads → C-allele  │
                           │ No cross-assignment │
                           └─────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                             COUNT GENERATION                                │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│    DNA Counts       │    │     RNA Counts      │    │  Count Statistics   │
│                     │    │                     │    │                     │
│ G-allele: 1000      │    │ G-allele: 500       │    │ • Barcode coverage  │
│ C-allele: 1000      │    │ C-allele: 2000      │    │ • Replicate quality │
│                     │    │                     │    │ • Count distribution│
│ Ratio_G = 500/1000  │    │ Ratio_C = 2000/1000 │    │                     │
│        = 0.5        │    │         = 2.0       │    │                     │
└─────────────────────┘    └─────────────────────┘    └─────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            VARIANT ANALYSIS                                 │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│  Individual Variant │    │  Master Variant     │    │   Correlation       │
│     Tables          │    │     Table           │    │    Analysis         │
│                     │    │                     │    │                     │
│ Per replicate:      │    │ Combined replicates:│    │ • Replicate r ≥ 0.8│
│ • REF counts        │    │ • Aggregated counts │    │ • Effect consistency│
│ • ALT counts        │    │ • Quality filtering │    │ • Statistical power │
│ • log2FC per rep    │    │ • Final log2FC      │    │                     │
└─────────────────────┘    └─────────────────────┘    └─────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          REGULATORY SNP CALLING                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│    Effect Size      │    │  Statistical Test   │    │  Final Results      │
│   Calculation       │    │                     │    │                     │
│                     │    │ • t-test/Mann-Whitney│    │ • Regulatory SNPs   │
│ log2FC = log2(2.0/0.5)   │ • Multiple testing  │    │ • Effect directions │
│        = 2.0        │    │ • FDR correction    │    │ • Confidence levels │
│                     │    │                     │    │                     │
│ 4-fold increase     │    │ p < 0.05, FDR < 0.1│    │ chr1:123 G→C        │
│ in C-allele         │    │                     │    │ Activating, 4x ↑    │
└─────────────────────┘    └─────────────────────┘    └─────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              OUTPUT FILES                                   │
│                                                                             │
│ • variant_table.tsv    - All variant measurements                          │
│ • correlation_table.tsv - Replicate correlations                           │
│ • qc_report.html       - Quality control metrics                           │
│ • regulatory_snps.tsv  - Significant regulatory variants                    │
└─────────────────────────────────────────────────────────────────────────────┘

KEY SUCCESS FACTORS:
═══════════════════════
1. EXACT MODE MAPPING → Eliminates cross-contamination (17% error prevention)
2. SUFFICIENT REPLICATES → ≥3 biological replicates for statistical power  
3. COUNT THRESHOLDS → min_dna_counts ≥10, min_rna_counts ≥5
4. QUALITY MONITORING → Replicate correlation r ≥ 0.8
5. EFFECT SIZE LIMITS → |log2FC| ≥ 0.3 for detection confidence
```