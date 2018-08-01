<p align="center"> 
<img src="logo/logo.png">
</p>

EnrichM is a set of tools for comparative genomics large sets of metagenome assembled genomes (MAGs). The current functionality is as follows:

1. A basic pipeline for annotating population genomes.
2. A pipeline to identify genes or metabolic pathways that are enriched within and between user-defined groups  
3. To construct metabolic networks from annotated population genomes. 

EnrichM is under active development, so use at your own risk. Currently, there are the following sub-commands:

```
Annotation 
	annotate		-> Basic annotation of bins and contigs.

Enrichment analysis
	classify		-> Determine what KEGG modules a genome encodes.
	enrichment		-> Generate an enrichment matrix from modules produced by 
					   annotate.

Network analysis
	pathway			-> Generate a metabolic network from specific KEGG module or 
					   compounds.
	explore			-> Explore a metabolic network from a given compound.
```

# Installation

## Dependencies 

EnrichM has the following non-python dependencies:
* [fxtract](https://github.com/ctSkennerton/fxtract) >= 2.3
* [hmmer](http://hmmer.org/) >= 3.1b
* [seqmagick](https://fhcrc.github.io/seqmagick/) >= 0.6.1
* [diamond](https://github.com/bbuchfink/diamond) >= 0.9.6
* [prodigal](http://prodigal.ornl.gov/) >= 2.6.3

## PyPi 
```
sudo pip install enrichm
```

# Setup
Before running enrichm, you'll need to download the back-end databases. This can be done using a command in enrichm. You may need super user privelages:
```
enrichm data
```
This should take approximately 15 minutes. To check for updates and install updates, simply run the same command. 

# Subcommands

## annotate
Annotate is a function that allows you to annotate your population genomes with [KO](http://www.kegg.jp/kegg/ko.html), [PFAM](http://pfam.xfam.org/), [TIGRFAM](http://www.jcvi.org/cgi-bin/tigrfams/index.cgi). The result will be a .gff file for each genome, and a frequency matrix for each annotation type where the rows are annotation IDs and the columns are genomes. 

See the [annotate help page](https://github.com/geronimp/enrichM/wiki/annotate) for more


## classify
Classify quickly reads in KO annotations in the form of a matrix (KO IDs as rows, genomes as columns) and determines which [KEGG modules](http://www.kegg.jp/kegg/module.html) are complete. Annotation matrices can be generated using the annotate function. 

See the [classify help page](https://github.com/geronimp/enrichM/wiki/classify) for more


## enrichment
Enrichment will read in KO or PFAM annotations in the form of a matrix (IDs as rows, genomes as columns) and a metadata file that separates genomes into groups to compare, and will run some basic stats to determine the enrichment of modules or pfam clans between and within the groups. 

See the [enrichment help page](https://github.com/geronimp/enrichM/wiki/enrichment) for more


## pathway
Pathway reads in a KO matrix and generates a Cytoscape-readable metabolic network and metadata file. Only reactions that are possible given the KOs present in the input matrix are shown, and the modules and reactions that are included in the output can be customized.

See the [pathway help page](https://github.com/geronimp/enrichM/wiki/pathway) for more


## explore
Explore is similar to pathway, but rather than generating a specified pathway it will start from a given query compound ID, and explore the possible reactions that use that compound given the enzymes present in the input KO matrix.

See the [explore help page](https://github.com/geronimp/enrichM/wiki/explore) for more

# Contact
If you have any feedback about EnrichM, drop an email to the [SupportM](https://groups.google.com/forum/?hl=en#!forum/supportm) public help forum. Software by [Joel A. Boyd](https://ecogenomic.org/personnel/mr-joel-boyd) (@geronimp) at the Australian Centre for Ecogenomics.

# License
EnrichM is licensed under the GNU GPL v3+. See LICENSE.txt for further details. 