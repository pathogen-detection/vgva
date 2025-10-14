# VGVA – Viral Genome Variation Analyzer

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Cite](https://img.shields.io/badge/cite-Gong%20et%20al.%202025-red)](https://doi.org/10.XXXX/XXXXXX)

---

## Overview

**VGVA** is a lightweight python package for the **whole-genome variation analysis** of DNA/RNA viruses.  
The virus whole-genome alignment matrix was segmented into consecutive sub-matrices (window=100bp, step=25bp) by sliding from left to right. 
For each sub-matrix, 
the python scripts will calculate the 25th percentile (Q1), median, and 75th percentile (Q3), upper bound, lower bound, and outliers of pairwise sequence similarities.

Results can be exported as tables or publication-ready SVG figures together with a user-supplied genome structure table.

&gt; If you use VGVA in your research :  
- Please cite our paper: Gong B, Xiang L, Li J, *et al.* *Novel Real-Time Quantitative RT-PCR for Detection of PRRSV-1 and PRRSV-2 Strains Circulating in China: A Study Based on Whole-Genome Characteristics and Epidemiological Features*  
- Obey the **GNU General Public License v3 (GPLv3)**.

---

## Installation

```bash
pip install vgva
```

- or clone and install in development mode:

```bash
git clone https://github.com/pathogen-detection/vgva.git
cd vgva
pip install -e . 
```
- or download the wheel file and install:
```bash
pip install dist/vgva-0.0.1-py3-none-any.whl
```

## Quick Start

```python
import os
import pandas as pd
from vgva import GenomeVariation
from vgva import set_text_editable

set_text_editable(True)  # keep SVG text selectable
if __name__ == "__main__":
    # 1. load alignment
    gv = GenomeVariation("test_data/EU02_mafft.fasta")
    
    # 2. compute variation metrics
    gv.get_genome_variation(
        window=100,  # window size (bp)
        step=25,  # slide step (bp)
        outliers_type="lower",  # report lower outliers
        gap_symbol="-",  # gap character
        is_include_gaps=True,  # gaps participate in similarity
        threads=os.cpu_count() * 2
    )
    
    # 3. load genome structure
    genome_structure = pd.read_excel("test_data/structure.xlsx", sheet_name="PRRSV-1")
    
    # 4. plot
    gv.plot_genome_variation(
        None,  # auto colour-map
        genome_structure_table=genome_structure,
        target_x_start=95,
        target_x_end=650,
        xaxis_step=200,
        yaxis_limit=(-0.05, 1.05),
        figures_height=(0.65, 0.2),  # (similarity, structure)
        figures_bottom=(0.31, 0.02),  # bottom offset in canvas
        outdir="test_data/",
        is_show=True
    )

```

## Genome-Structure Table
Excel (or CSV) file with mandatory columns:
### Genome-Structure Table (Excel, CSV, TSV, AND THE LIKE)

Required columns (order does **not** matter):

| Column        | Type  | Description                                   |
|---------------|-------|-----------------------------------------------|
| `gene`        | str   | Feature name (e.g. *ORF5*, *5'UTR*)           |
| `start_pos`   | int   | 1-based start position on the genome          |
| `end_pos`     | int   | 1-based end position (≥ start_pos)            |
| `orientation` | int   | 1 = forward / sense; -1 = reverse / antisense |
| `y_offset`    | float | gene name y offset (0.0, 1.0)                 |
| `x_offset`    | float | gene name x offset (0.0, 1.0)                 |

Example snippet(PRRSV-2):

| gene  | start_pos | end_pos | orientation | y_offset | x_offset |
|-------|-----------|---------|-------------|----------|----------|
| 5'UTR | 1         | 190     | 1           | 0.6      | -0.3     |
| ORF1a | 191       | 7804    | 1           | 0.6      | 0.5      |
| ORF1b | 7786      | 12174   | -1          | 0.8      | 0.5      |
| ORF2a | 12176     | 12946   | 1           | 0.6      | 0        |
| ORF2b | 12181     | 12402   | -1          | 0.8      | 0        |
| ORF3  | 12799     | 13563   | -1          | 0.8      | 0        |
| ORF4  | 13344     | 13880   | 1           | 0.6      | 0        |
| ORF5  | 13891     | 14493   | -1          | 0.8      | 0        |
| ORF5a | 13881     | 14036   | 1           | 0.6      | 0        |
| ORF6  | 14478     | 15002   | 1           | 0.6      | 0        |
| ORF7  | 14992     | 15363   | -1          | 0.8      | 0        |
| 3'UTR | 15364     | 15514   | 1           | 0.6      | 0        |
| PolyA | 15515     | 15559   | -1          | 0.8      | 0        |
| 1α    | 191       | 730     | 1           | 0.2      | 0.5      |
| 1β    | 731       | 1336    | 1           | 0.2      | 0.5      |
| 2     | 1337      | 5029    | 1           | 0.2      | 0.5      |
| 3     | 5030      | 5719    | 1           | 0.2      | 0.5      |
| 4     | 5720      | 6331    | 1           | 0.2      | 0.5      |
| 5     | 6332      | 6841    | 1           | 0.2      | 0.5      |
| 6     | 6842      | 6889    | 1           | 0.6      | 0        |
| 7α    | 6890      | 7336    | 1           | 0.2      | 0.2      |
| 7β    | 7337      | 7666    | 1           | 0.2      | 0.2      |
| 8     | 7667      | 7801    | 1           | 0.2      | 0        |
| 9     | 7786      | 9720    | -1          | 0.3      | 0.5      |
| 10    | 9721      | 11043   | -1          | 0.3      | 0.4      |
| 11    | 11044     | 11712   | -1          | 0.3      | 0.4      |
| 12    | 11713     | 12174   | -1          | 0.3      | 0.2      |

## Requirements
- Python 3.8+
- numpy>=1.24.4
- pandas>=2.0.3
- matplotlib>=3.6.2
- biopython>=1.83
