# Step 14: external metadata gate (preliminary source registry)

This is a source-screening plan, not an external model result. The initial registry is
[`metadata/v01_external_catalog/candidates_preliminary.tsv`](../metadata/v01_external_catalog/candidates_preliminary.tsv).
All values marked “待核验/待计算” remain unknown. Do not equate cells, GEM
batches, donor or biological replicates, and effective perturbation tasks.

## Verified source records

- **Jiang et al.** Primary [paper](https://doi.org/10.1038/s41556-025-01622-z)
  reports >1,500 perturbations across six cell lines and five signalling contexts.
  The [GEO record GSE281048](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281048)
  names A549, MCF7, HT29, HAP1, BxPC3 and K562; IFNB, IFNG, TGFB,
  TNFA and INS; 44–61 selected regulators per pathway, three independent guides
  per gene, 14 non-targeting controls per pathway library, and four pooled GSM
  records marked by replicate. “>1,500 perturbations” is **not** a confirmed
  count of distinct target genes shared with our 1,264-task panel. GEO says
  its Seurat supplementary files contain raw counts and cell metadata, but
  donor/batch/guide usability needs direct inspection.
  The authors' [Zenodo v2.1 record](https://zenodo.org/records/14518762)
  lists a 1.5 kB README, a 324.1 MB DE ZIP and five much larger pathway RDS
  files; the displayed collection size is about 20.5 GB.
- **X-Atlas/Orion.** The [preprint](https://doi.org/10.1101/2025.06.11.659105)
  reports two screens (HCT116 and HEK293T) with about eight million cells.
  The [author dataset card](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/blob/main/README.md)
  says 18,903 protein-coding genes were targeted, median ≥140 cells per
  perturbation, on-target knockdown medians of 75.4% and 51.5% respectively,
  and `sample` means GEM batch. It lists `gene_expression` as sparse raw
  counts and `gene_token_id` as the link to Ensembl IDs.
  The [author Figshare+ v3 record](https://doi.org/10.25452/figshare.plus.29190726)
  lists 521.24 GB of **processed** data, dual-guide-filtered H5ADs, guide
  library CSV and guide-call CSV files. The author
  [Hugging Face tree](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/tree/main)
  shows about 126 GB of Parquet files and an
  [882 kB gene metadata file](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/tree/main/metadata).
  The author card and Figshare+ page display CC BY-NC-SA 4.0.
  Eight million cells do **not** imply independent biological replicates
  or sufficient GEM-matched non-targeting cells.

## Gate before any full single-cell file

1. Save source metadata, sample/condition/replicate tables, perturbation and
   guide labels, control definitions, file sizes and hashes. Restrict early
   download to small metadata. The Jiang DE ZIP must not be inspected for
   effect sizes while choosing a target or risk rule.
2. Match external target symbols to our 1,264 target labels and external
   Ensembl IDs to 6,642 expression gene IDs. Calculate counts **per cell
   line × stimulus × batch × perturbation × guide**, not just overall.
3. Calculate non-targeting cells in each matching batch and confirm guide,
   replicate and knockdown definitions. Inspect confounding of perturbation
   with batch, donor and stimulus. Treat GEM as a technical grouping unless
   independent biological repeats are documented.
4. Record overlap with the four existing raw datasets, including K562 in
   Jiang, and the guide design shared between X-Atlas and Replogle. Decide
   on one eligible external *target* context before examining its response.
5. Freeze prediction, calibration, risk rule, failure definition, zero-response
   baseline, metrics, code commit and SHA256 before using any external target
   expression, pseudobulk or DE response. Record the freeze timestamp.
   Prior HepG2 analysis is retrospective, not independent validation.

## Next metadata-only source: Jiang guide design

The article [Supplementary Tables workbook](https://www.nature.com/articles/s41556-025-01622-z)
contains Supplementary Tables 1–5. The manuscript identifies Table 1 as the
pathway sgRNA and non-targeting control design. Download the workbook from the
official article's “Supplementary information” section; record the final URL,
byte size and SHA256. Inspect **Table 1 only** to derive the designed guide-to-
target/control map and potential overlap with the project's 1,264 labels.
Keep other sheets unopened until their content has been classified as safe
metadata or response-derived outcomes. Guide design does not establish observed
per-guide cell support, GEM-matched controls, sample replication or expression
feature overlap. Those still require separate cell metadata and raw count
checks. Do not inspect the Zenodo DE ZIP for cohort selection.

## Observed Jiang supplementary workbook structure (HPC metadata check)

The operator downloaded the official Supplementary Tables XLSX as
`41556_2025_1622_MOESM3_ESM.xlsx` (2,962,962 bytes). Its ZIP integrity
check passed. The workbook has `Table_of_Contents`, ST1a–ST1e guide
sheets (IFNG, IFNB, TGFB, TNFA, INS), `ST1f_all_genes_by_pathway`,
`ST2_cytokine_concentrations`, `ST3_Outliers_Ultima_Illumina`,
ST4a–ST4d MultiCCA sheets and ST5a–ST5e hclust sheets. Its SHA256
has not yet been pasted into the audit record. Before a model freeze,
inspect only ST1a–ST1f design labels and headers. Do not read ST3–ST5
analytical outcomes to select the validation cohort or risk rules.

## First operator action

Check the local code repository is clean and fetch only the small Zenodo
catalog JSON; no large RDS/H5AD/Parquet transfer:

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
mkdir -p metadata/v01_external_catalog
curl --fail --location --retry 2 --max-time 60 \
  "https://zenodo.org/api/records/14518762" \
  -o metadata/v01_external_catalog/jiang_zenodo_record.json
python3 -m json.tool metadata/v01_external_catalog/jiang_zenodo_record.json >/dev/null
sha256sum metadata/v01_external_catalog/jiang_zenodo_record.json
```

If HTTP or DNS access fails, preserve the curl error and inspect the network
policy before trying another data endpoint. Do not interpret download failure
as absence of public data. Review the JSON file listing/rights before any
larger download and add the subsequent actual project-computed coverage to
a new versioned result.
