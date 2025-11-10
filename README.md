# 🧬 Viractin Detection & Phylogenetic Analysis (Anvi'o 7.1)

This repository contains the complete workflow used in the **NCLDV Viractin TD/TP** practical,
where we identify *actin-related proteins (Viractins)* from metagenome-assembled genomes (MAGs) of large dsDNA viruses (NCLDVs)
and reconstruct their phylogenetic relationships.

---

## 📂 Project Overview

| Step | Description | Key Command |
|------|--------------|--------------|
| **1. Create MAG Databases** | Generated Anvi’o databases for each MAG using contigs | `anvi-gen-contigs-database -f data/<MAG>.fa -o myData/<MAG>-DB.db` |
| **2. Functional Annotation (PersoHMM)** | Added custom HMM profiles for actin-related genes | `anvi-run-hmms -c NCLDV-GENOMES.db -H hmms/PersoHMM/` |
| **3. Gene Calling & QC** | Evaluated genome completeness and redundancy | `anvi-display-contigs-stats myData/*-DB.db --report-as-text` |
| **4. Extract AA Sequences** | Exported amino acid sequences for all genes | `anvi-get-sequences-for-gene-calls --get-aa-sequences` |
| **5. HMM Search** | Detected actin-like hits with HMMER | `hmmscan --tblout results/hmmer/<MAG>.tbl -E 1e-12 hmms/PersoHMM/genes.hmm` |
| **6. Sequence Extraction** | Retrieved hits using seqkit | `seqkit grep -f <ids> <aa.fa>` |
| **7. Multiple Sequence Alignment** | MUSCLE / trimAl for cleaned alignment | `muscle -in actin_plus_refs.fa -out actin_plus_refs.muscle.fa` |
| **8. Phylogenetic Tree** | IQ-TREE maximum likelihood with LG+R7 | `iqtree -s actin_plus_refs.trimal.phy -m LG+R7 -alrt 1000 -bb 1000` |

---

## 🧠 Results Summary

- **20 MAGs analyzed**
- **HMM source:** PersoHMM (custom actin-like model)
- **Hits:** 1 Viractin sequence per MAG
- **Alignment:** 368 amino acids across all sequences
- **Model:** LG+R7
- **Bootstrap:** 1000 ultrafast + SH-aLRT 1000
- **Output tree:** `results/actin_plus_refs.trimal.phy.treefile`

---

## 🌳 Phylogenetic Tree

The final tree (`actin_plus_refs.trimal.phy.contree`) was visualized with **ete3** / **FigTree**.

<p align="center">
  <img src="results/actin_plus_refs_tree_support.png" width="600">
</p>

---

## 🧪 Tools & Environment

| Tool | Version | Source |
|------|----------|--------|
| Anvi’o | 7.1 | bioconda |
| HMMER | 3.3.2 | bioconda |
| MUSCLE | 3.8.1551 | bioconda |
| trimAl | 1.4.1 | bioconda |
| IQ-TREE | 1.6.12 | bioconda |
| BMGE | 1.12 | bioconda |
| seqkit | 2.3.1 | bioconda |

---

## 🧰 Python Visualization Code (Optional)

This script visualizes your final tree with branch support using **ete3**:

```python
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

tree = Tree("results/actin_plus_refs.trimal.phy.contree")

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_support = True
ts.scale =  50
ts.title.add_face(TextFace("Viractin Phylogeny (LG+R7 Model, 20 Sequences)", fsize=14), column=0)

nstyle = NodeStyle()
nstyle["fgcolor"] = "#444"
nstyle["size"] = 0
for n in tree.traverse():
    n.set_style(nstyle)

tree.render("results/actin_plus_refs_tree_support.png", tree_style=ts)

