Detecting **novel bacteria or viruses**—those not present in existing reference databases—is a complex but increasingly tractable challenge thanks to advances in **NGS**, **assembly algorithms**, and **ML/DL models**. 

I created  a **roadmap** that combines metagenomics, de novo analysis, and advanced machine learning to maximize chances of detecting and characterizing unknown pathogens.

---

## 🧭 **Roadmap for Detecting Novel Bacteria or Viruses from NGS Data**

### **1. Sample Preparation**

* **Source**: Clinical or pharmaceutical (tissue, blood, swab, product, environmental).
* **Preservation**: Ensure RNA/DNA integrity (e.g., RNAlater, cold chain).
* **Enrichment (optional)**: Use microbial enrichment (e.g., DNase treatment to remove host DNA).

---

### **2. Preprocessing**

| Step                     | Tools                                                 |
| ------------------------ | ----------------------------------------------------- |
| **Quality Control**      | `FastQC`, `MultiQC`, `Trimmomatic`                    |
| **Host Read Removal**    | `Bowtie2` against human genome (`hg38`) or `BMTagger` |
| **rRNA Depletion Check** | `SortMeRNA`, `Barrnap`                                |

---

### **3. De Novo Assembly**

To detect **novel genomes**, reference-free assembly is critical.

* **Tools**: `metaSPAdes`, `MEGAHIT`, `SPAdes`, `Flye` (long-read)
* **Assembly QC**: `QUAST`, `CheckM`, `MetaQUAST`

---

### **4. Contig Binning**

Group contigs into species-like genome bins (MAGs = metagenome-assembled genomes).

* **Tools**: `MetaBAT2`, `MaxBin2`, `CONCOCT`, `VAMB` (uses VAE + clustering)
* **Quality Control**: `CheckM`, `BUSCO`

---

### **5. Taxonomic Assignment**

**Compare against reference databases**, but focus on **unclassified/low-identity hits**.

* **Tools**:

  * `Kraken2`, `Kaiju`, `Centrifuge`, `MetaPhlAn` (for known)
  * `CAT/BAT`, `MMseqs2`, `Mash`, `Sourmash` for similarity-based novelty screening

* **Strategy**:

  * Identify contigs with **no high-identity match** to RefSeq or GTDB
  * Flag <90% ANI (Average Nucleotide Identity) as **potentially novel**

---

### **6. Gene Prediction & Functional Annotation**

Even novel genomes have conserved domains:

* **Gene prediction**: `Prodigal`, `MetaGeneMark`
* **Functional annotation**: `eggNOG-mapper`, `PfamScan`, `InterProScan`, `DRAM`
* **Virus-specific tools**: `VirSorter2`, `VIBRANT` to identify viral elements

---

### **7. ML/DL Novelty Detection Pipeline**

**This is key to discovering *truly novel* organisms**, especially viruses or divergent bacteria:

#### **a. Embedding Sequences**

* Use **transformers** like `DNABERT`, `ESM`, `ProtT5` to embed sequences without reference.

#### **b. Anomaly Detection**

* Apply **unsupervised learning** to detect outliers from known distributions:

  * `Isolation Forest`, `Autoencoders`, `DeepSVDD`
  * Train on known microbes → Flag unknown-like behavior (e.g., coding potential, GC content, k-mer composition)

#### **c. Clustering & Visualization**

* Tools: `t-SNE`, `UMAP`, `HDBSCAN`
* Cluster contigs/MAGs based on sequence embeddings or functional profiles.

---

### **8. Viral/Bacterial Specific Identification**

* **Bacteria**: Use GTDB-Tk for phylogenetic placement and determine ANI for novelty.
* **Viruses**:

  * `CheckV`, `VirFinder`, `DeepVirFinder`: Detect and assess viral completeness.
  * `GraphBin`, `vConTACT2`: Classify viral contigs and infer networks of relatedness.

---

### **9. Validation & Downstream Investigation**

* **Confirm novelty** via phylogenetic analysis: `FastTree`, `IQ-TREE`, `PhyloPhlAn`
* Submit to **NCBI GenBank**, **ENA**, or **IMG/M** for comparison
* Perform **epidemiological screening** on more samples

---

### **10. Regulatory & Biosafety Check (U.S. Context)**

* Flag any novel pathogen under **select agent regulations** (CDC, USDA)
* Follow **Biosafety Level (BSL)** guidance for handling unknowns
* Consider IRB and **Dual Use Research of Concern (DURC)** classification

---


