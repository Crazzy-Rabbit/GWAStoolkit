# 🧬 GWAStoolkit
**A unified, high-performance C++ toolkit for processing GWAS summary statistics.**

GWAStoolkit integrates multiple commonly needed GWAS operations into a single, efficient command-line tool. \
It provides consistent interfaces, shared parameters across subcommands, and fast performance for very large datasets (dbSNP > 15GB, millions of SNPs). 

## Features
### 1️⃣ rsidImpu — Fast rsID annotation using dbSNP

- Supports TSV / CSV / gzipped input
- Allele-aware matching with:
    - A1/A2 swapping
    - Strand complement (A↔T, C↔G)
- dbSNP / bim supported
- Optional output formats (COJO, POPCORN, MR-MEGA, etc.)
- Automatic QC: MAF, beta, se, p, freq, N
- Remove duplicate SNPs by smallest P-value

### 2️⃣ convert — Convert between GWAS formats

Convert any GWAS summary file to formats required by:
- **GWAS standard**
- **GCTA-COJO**
- **POPCORN**
- **MR-MEGA**

### 3️⃣ or2beta — Convert OR → beta + SE

- Converts OR to log-odds beta
- Computes SE from OR, or from P-value if needed
- Full QC support

### 4️⃣ computeNeff — Compute effective sample size (binary traits)

Used for case/control GWAS.
- Full QC support
- Supports two modes:
   - **Mode 1 — Fixed case/control numbers**
   - **Mode 2 — Per-SNP case/control columns**
```
Neff = 4 * case * control / (case + control)
```
After computing Neff, all SNPs are standardized:
```
z   = beta / se
se  = 1 / sqrt(2 * p * (1 - p) * (Neff + z^2))
beta= z * se
```
## ⚙️ Installation
### Dependencies
- g++ (support C++11)
- zlib
### Build
```
git clone https://github.com/Crazzy-Rabbit/GWAStoolkit.git
cd GWAStoolkit
make
```
You will get:
```
./GWAStoolkit
```
## 🚀 Quick Start
### List all commands
```
GWAStoolkit --help
```
### 🔧 1. rsidImpu Example

Basic:
```
GWAStoolkit rsidImpu \
  --gwas-summary gwas.txt.gz \
  --dbsnp dbsnp.txt.gz \
  --out gwas.rsid.txt.gz \
  --dbchr CHR --dbpos POS --dbA1 REF --dbA2 ALT --dbrsid RSID
```

With specific output format:
```
--format cojo \
--freq Freq --beta Beta --se SE --n N --pval P
```
### 🔧 2. convert Example
```
GWAStoolkit convert \
  --gwas-summary gwas.txt \
  --out gwas.cojo.txt \
  --format cojo \
  --SNP SNP --A1 A1 --A2 A2 --pval P \
  --freq FRQ --beta BETA --se SE --n N
```
### 🔧 3. or2beta Example
```
GWAStoolkit or2beta \
  --gwas-summary gwas.txt \
  --out gwas.beta.txt \
  --or OR \
  --SNP SNP --A1 A1 --A2 A2 --pval P
```
### 🔧 4. computeNeff Example

Mode 1 — global case/control:
```
GWAStoolkit computeNeff \
  --gwas-summary gwas.txt \
  --out gwas.neff.txt \
  --case 20000 --control 30000
```

Mode 2 — per-SNP:
```
GWAStoolkit computeNeff \
  --gwas-summary gwas.txt \
  --out gwas.neff.txt \
  --case-col CASE \
  --control-col CTRL
  ```
## 📦 Unified Argument System

All subcommands share:
| Parameter        | Description                                  | Default          |
|------------------|----------------------------------------------|------------------|
| `--gwas-summary` | Input GWAS (txt/tsv/csv/gz)                  | required         |
| `--out`          | Output file (txt/gz supported)               | required         |
| `--format`       | gwas/cojo/popcorn/mrmega                     | gwas             |
| `--chr` `--pos` `--A1` `--A2` `--pval` | Column names           | CHR/POS/A1/A2/p  |
| `--freq` `--beta` `--se` `--n` | Effect model columns            | freq/b/se/N      |
| `--maf`          | MAF threshold                                | 0.01             |
| `--remove-dup-snp` | Drop duplicated SNP (keep smallest P)      | off              |
| `--threads`      | Multi-threading                              | 1                |
| `--log FILE`     | Write log file                               | none             |

Additional command-specific parameters:

| Command      | Extra Required Parameters                          |
|--------------|---------------------------------------------------|
| rsidImpu     | `--dbsnp --dbchr --dbpos --dbA1 --dbA2 --dbrsid` |
| convert      | `--SNP`                                           |
| or2beta      | `--SNP --or`                                     |
| computeNeff  | `(--case & --control) or (--case-col & --control-col)` |

## 🧪 Example Output (COJO Format `--format cojo`)

```
SNP       A1  A2  freq   b    se      p       N
rs1000    A   G   0.37   0.145   0.035   1e-5    50000
rs2000    T   C   0.42  -0.080   0.025   2e-3    50000
...
```