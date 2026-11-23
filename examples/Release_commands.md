

## **1. First pass — November 11 release**

The first database should represent the original release only:

| Component            | First-pass input                    |
| -------------------- | ----------------------------------- |
| Cultured baseline    | Hungate                             |
| Secondary baseline   | Rumen_Isolates                      |
| Marker cohorts       | QUB, UoG, BSI                       |
| Genome-only isolates | AAFC (34) + IAS_CAAS (878)          |
| Genome results       | 912 existing genome results         |
| Marker reviews       | The original reviewed decision file |
| GTDB                 | r232                                |


### **Step 1 — Build a fresh database**

From your November Release working directory, using the cleaned project config:

```
branchmanager filing-cabinet \
  --config Rumen_Gateway_Assistant.json \
  --db project.sqlite \
  --build-tree \
  --threads 10 \
  -o 01_filing_cabinet
```

This creates the fresh database and loads the configured Hungate and Rumen_Isolates baselines.

Do this in a **new run directory**, not the existing November 19 database.

### **Step 2 — Register the genome-only isolates**

Alternatively, a genome-results TSV containing `sequence_id` and `partner_id`
can register all listed isolates directly, using one candidate dataset:

```
branchmanager import-partner-workbook \
  --db project.sqlite \
  --input Genome_Results.tsv \
  --dataset Genome_Only \
  --genome-only \
  -o 02_genome_registration
```

Registration imports identities only, not genome assemblies, QC or taxonomy.
Existing isolates assigned to a different dataset are reported as conflicts and
are not reassigned. Import the genome results separately in Step 3.

The first-pass workbook contains both AAFC and IAS_CAAS. As before, it needs to be split into two partner-specific workbooks so each has the correct dataset designation.

Then run:

```
branchmanager import-partner-workbook \
  --db project.sqlite \
  --input AAFC_Genome_Only_Registration.xlsx \
  --sheet "Bulk Upload Template" \
  --dataset AAFC_01 \
  --genome-only \
  -o 02a_AAFC_registration
```

```
branchmanager import-partner-workbook \
  --db project.sqlite \
  --input IAS_CAAS_Genome_Only_Registration.xlsx \
  --sheet "Bulk Upload Template" \
  --dataset IAS_CAAS_01 \
  --genome-only \
  -o 02b_IAS_CAAS_registration
```

### **Step 3 — Import the existing genome results**

```
branchmanager import-genome-results \
  --db project.sqlite \
  --results-tsv AAFC_IAS_CAAS_Genome_Results.tsv \
  --out 03_genome_results_import
```

This imports the supplied QC and GTDB taxonomy results. It does not rerun genome analysis.

### **Step 4 — Prepare the three marker cohorts**

Use the November 11/first-release catalogue and manifest:

```
branchmanager assistant \
  --config Rumen_Gateway_Assistant.json \
  --db project.sqlite \
  --partner-metadata Partner_Isolate_Catalogue_November.csv \
  --cohort-manifest Partner_Marker_Cohort_Manifest.csv \
  --marker-review cohort_marker_review_decisions_clean_20261105.tsv \
  --cohort-prepare-only \
  --threads 10 \
  --out 04_cohort_preparation
```

This performs the preparation/QC stage and generates the cohort review materials. It does **not** publish the final joint assessment to the database yet.

### **Step 5 — Review and reconcile**

Before resuming:

- Confirm the expected QUB, UoG and BSI isolate IDs are present.
- Check the generated review template.
- Confirm that the historical review decisions match the generated sequence IDs.
- Investigate orphaned decisions and newly unreviewed sequences.
- Do not automatically approve new or unmatched sequences.

I would retain the review reconciliation step from the old script as a small separate utility/reporting step rather than reintroducing the entire Bash workflow.

### **Step 6 — Resume the prepared cohort**

Once the decisions are confirmed and reconciled:

```
branchmanager assistant \
  --config Rumen_Gateway_Assistant.json \
  --db project.sqlite \
  --partner-metadata Partner_Isolate_Catalogue_November.csv \
  --cohort-resume 04_cohort_preparation \
  --marker-review matched_marker_review_decisions.tsv \
  --threads 10 \
  --out 04_cohort_preparation
```

Important: The exact `matched_marker_review_decisions.tsv` filename depends on your reconciliation output. Do not use that filename unless the reconciliation step actually generated it. The review file must be verified against the freshly generated template before resume.

## **3. Continuation — November 20 release**

The November 20 inputs are cumulative rather than a replacement of the first release.

## **What changes from the first pass?**

**November 11**

1,788 catalogue records

QUB · UoG · BSI · AAFC · IAS_CAAS

**November 20**

1,849 catalogue records

\+ AU · UNA · UoM additions

Based on the uploaded files, the November 20 update contains:

| Partner   | Added records relative to Nov 11 |
| --------- | -------------------------------- |
| AU        | 39                               |
| UNA       | 2                                |
| UoM       | 20                               |
| **Total** | **61**                           |

The two new marker partners are **AU and UNA**. The additional 20 UoM records are genome-only additions in the updated genome results.

### **Step 7 — Continuation dry run**

Once the first-pass database is validated, use the cumulative November 20 catalogue and genome results:

```
branchmanager continue-project \
  --db project.sqlite \
  --partner-catalogue Partner_Isolate_Catalogue_November.csv \
  --genome-results Genome_Results.tsv \
  --partner-workbook AarhusUniversity_DATA_2026-11-18_1017.csv \
  --partner-workbook-dataset AU_01 \
  --dry-run \
  -o 08_continue_nov20_dryrun
```

The `--partner-workbook` here is the AU CSV. The continuation workflow accepts CSV metadata; the sheet option is unnecessary for it.

### **Step 8 — Continuation apply**

Only after the dry-run reports are satisfactory:

```
branchmanager continue-project \
  --db project.sqlite \
  --partner-catalogue Partner_Isolate_Catalogue_November.csv \
  --genome-results Genome_Results.tsv \
  --partner-workbook AarhusUniversity_DATA_2026-11-18_1017.csv \
  --partner-workbook-dataset AU_01 \
  -o 08_continue_nov20
```

This step updates the cumulative catalogue and genome evidence. **It does not itself process the new AU and UNA marker sequences.**

### **Step 9 — Process the two new marker partners**

For the additional marker work, I would create a new manifest containing **only AU and UNA** from the November 20 five-partner manifest. I would not rerun QUB, UoG and BSI unnecessarily.

Then prepare the new marker cohort against the updated database, using the November 20 cumulative catalogue:

```
branchmanager assistant \
  --config Rumen_Gateway_Assistant.json \
  --db project.sqlite \
  --partner-metadata Partner_Isolate_Catalogue_November.csv \
  --cohort-manifest AU_UNA_Marker_Manifest.csv \
  --cohort-prepare-only \
  --threads 10 \
  --out 09_AU_UNA_preparation
```

Review the new marker sequences and QC decisions, then resume using the generated review file:

```
branchmanager assistant \
  --config Rumen_Gateway_Assistant.json \
  --db project.sqlite \
  --partner-metadata Partner_Isolate_Catalogue_November.csv \
  --cohort-resume 09_AU_UNA_preparation \
  --marker-review AU_UNA_review_decisions.tsv \
  --threads 10 \
  --out 09_AU_UNA_assessment
```

The review filename above is illustrative — use the actual generated and completed file.

## **4. Do we need code changes?**

### **For the first pass: No workflow code changes are inherently required**

The CLI already supports:

- Configured multiple baselines
- Fresh database construction
- Partner workbook registration
- Genome-results import
- Cohort manifest processing
- Preparation-only runs
- Cohort resumption

We can replace the large Bash script with a sequence of CLI commands and a small workbook-splitting/reconciliation utility where needed.

### **For the November 20 continuation: Yes, I would fix the continuation implementation first**

There are two important limitations in the current release I inspected:

1. `continue-project`**&#x20;handles catalogue and genome-result updates, but not marker-cohort processing.** That is why AU and UNA need a separate Assistant cohort workflow.
2. The latest BranchManager archive I inspected does **not contain the&#x20;**`continuation.py`**&#x20;module** in its source package. This is consistent with the packaging issue you previously encountered. I would not trust a continuation run from that release until this is corrected and tested.

So I would not run the November 20 continuation on your newly generated first-pass database until we have repaired and tested the installed package.

## **5. The workflow I recommend**

### **Release 1 — November 11**

Build the foundation

JSON project profile → Filing Cabinet (Hungate + Secondary) → AAFC/IAS_CAAS registration → genome-results import → QUB/UoG/BSI cohort preparation → review → resume.

### **Release 2 — November 20**

Continue the project cumulatively

Continuation dry-run → catalogue and genome update → AU/UNA marker preparation → new QC review → cohort resume → final validation.

This is preferable to the old Bash script because each operation is explicit, independently auditable, and easier to rerun if a particular stage fails. The JSON becomes a stable project profile rather than a release-specific script configuration.

**My recommendation is to do one more step before you start running this:** I can prepare the cleaned, partner-independent JSON and a concise command-runbook with the exact commands, the small AAFC/IAS workbook split utility, the AU/UNA manifest extraction, and validation checkpoints — all matched to the actual November 11 and November 20 files. That would give you a reproducible workflow without another monolithic Bash script.