<p align="center">
<img src="logo/logo.png">
</p>

[![Publish to PyPI](https://github.com/geronimp/enrichM/actions/workflows/publish.yml/badge.svg)](https://github.com/geronimp/enrichM/actions/workflows/publish.yml)

EnrichM is a set of comparative genomics tools for large sets of metagenome assembled genomes (MAGs). The current functionality includes:

1. A basic annotation pipeline for MAGs.
2. A pipeline to determine the metabolic pathways that are encoded by MAGs, using KEGG modules as a reference (although custom pathways can be specified).
3. A pipeline to identify genes or metabolic pathways that are enriched within and between user-defined groups of genomes (groups can be genomes that are related functionally, phylogenetically, recovered from different environments, etc).
4. Construct random forest machine learning models from the functional composition of MAGs, metagenomes or transcriptomes.
5. Apply random forest models to classify new MAGs or metagenomes.

EnrichM is under active development, so there is no guarantee that master is stable. It is recommended to install from a tagged release (see below).

# Installation
## Dependencies
EnrichM is written in Python 3 and requires >= 3.8. EnrichM requires the following non-Python dependencies:
* [hmmer](http://hmmer.org/) >= 3.4
* [diamond](https://github.com/bbuchfink/diamond) >= 2.0
* [prodigal](http://prodigal.ornl.gov/) >= 2.6.3
* [parallel](https://www.gnu.org/software/parallel/) >= 20180222
* [mmseqs2](https://github.com/soedinglab/MMseqs2) >= 13

## conda (recommended)
Clone the repository and create the conda environment:
```
git clone https://github.com/geronimp/enrichM.git
cd enrichM
conda env create -f environment.yml
conda activate enrichm
pip install .
```

## PyPI
```
pip install enrichm
```
Note: non-Python dependencies (hmmer, diamond, prodigal, parallel, mmseqs2) must be installed separately when using PyPI.

After installation, you'll need to download the back-end databases.

# Setup
## Loading EnrichM's database
The database contains Pfam-A HMMs, TIGRfam HMMs, dbCAN HMMs, and KoFamKOALA HMMs. By default it is installed in `~/enrichm_data`. Build it using:
```
enrichm data --create
```
To store the database in a custom location:
```
enrichm data --create --db_path /path/to/database/
```
To uninstall:
```
enrichm data --uninstall
```

## Using an existing database
If the database was built in a custom location, set the `ENRICHM_DB` environment variable so EnrichM can find it:

```
export ENRICHM_DB=/path/to/database/
```

Add this to your `.bashrc` or conda `activate.d` script to avoid setting it each session.

# Subcommands
## annotate
Annotate population genomes with [KO HMMs](http://www.kegg.jp/kegg/ko.html), [Pfam](http://pfam.xfam.org/), [TIGRfam](http://www.jcvi.org/cgi-bin/tigrfams/index.cgi), and CAZymes using [dbCAN](https://bcb.unl.edu/dbCAN2). The result is a GFF file for each genome and a frequency matrix for each annotation type (annotation IDs as rows, genomes as columns).

## classify
Reads KO annotations in the form of a matrix and determines which [KEGG modules](http://www.kegg.jp/kegg/module.html) are complete. Annotation matrices can be generated using `annotate`.

## enrichment
Enrichment reads an annotation matrix (IDs as rows, genomes as columns) and a metadata file separating genomes into groups, and runs statistical tests (Mann-Whitney U, Fisher's exact, Kruskal-Wallis) to identify enriched annotations between groups. Outputs include effect sizes, fold changes, and FDR-corrected p-values. Additional features include:

- **Synteny analysis**: identifies conserved gene blocks (operons) among enriched genes using intergenic distance thresholds
- **Mobile element proximity**: flags enriched genes located near transposases or insertion sequences
- Accepts output from `annotate`, or external tools including DRAM, eggNOG-mapper

## generate
Trains a random forest classifier or regressor from an annotation matrix and a metadata file of labels. Performs automated hyperparameter tuning via RandomizedSearchCV (with optional GridSearchCV refinement). Outputs the trained model, feature importances, and accuracy summary.

## predict
Applies a trained model (from `generate`) to a new annotation matrix and outputs per-sample predictions and class probabilities.

# Contact
If you have any feedback about EnrichM, drop an email to the [SupportM](https://groups.google.com/forum/?hl=en#!forum/supportm) public help forum. Software by [Joel A. Boyd](https://ecogenomic.org/personnel/mr-joel-boyd) (@geronimp) at the Australian Centre for Ecogenomics (ACE).

# License
EnrichM is licensed under the GNU GPL v3+. See LICENSE.txt for further details.

# Contributing 
I want EnrichM to be as useful as possible, so please feel free to leave feature requests and bug reports.

# Citation
If you find EnrichM useful and use it in your work, please cite it as follows:
```
Comparative genomics using EnrichM. Joel A Boyd, Ben J Woodcroft, Gene W Tyson. In preparation.
```
