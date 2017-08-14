![logo](logo/logo.png)

enrichM has three main purposes
1. To provide a very basic pipeline for annotating population genomes and assemblies. 
2. To compare user-defined groups of genomes and identify genes or metabolic pathways that are enriched within and between groups. 
3. To construct metabolic networks from annotated population genomes. 

enrichM is currently in a state of flux and currenty does not have unit tests so while you're welcome to use it, do so at your own risk.

enrichM currently has the following sub-commands:

```
  Genome annotation
    annotate    -> Annotate contigs with KO, PFAM, and TIGRfam IDs

  Enrichment analysis
    classify    -> Interpret KO annotations and report completeness of KEGG modules that a genome encodes
    enrichment  -> Compare groups of genomes and determine what groups of genes or functions are enriched in one set vs another.
    
  Network analysis
    pathway     -> Generate a Cytoscape-readable metabolic network from specific KEGG module or compounds.
    explore     -> Explore a metabolic network from a given compound            
```

# annotate
Annotate is a function that allows you to annotate your population genomes with [KO](http://www.kegg.jp/kegg/ko.html), [PFAM](http://pfam.xfam.org/), [TIGRFAM](http://www.jcvi.org/cgi-bin/tigrfams/index.cgi). The result will be a .gff file for each genome, and a frequency matrix for each annotation type where the rows are annotation IDs and the columns are genomes. 

See the [annotate help page](https://github.com/geronimp/enrichm/wiki) for more


# classify
Classify quickly reads in KO annotations in the form of a matrix (KO IDs as rows, genomes as columns) and determines which [KEGG modules](http://www.kegg.jp/kegg/module.html) are complete. Annotation matrices can be generated using the annotate function. 

See the [classify help page](https://github.com/geronimp/enrichm/wiki) for more


# enrichment
Enrichment will read in KO or PFAM annotations in the form of a matrix (IDs as rows, genomes as columns) and a metadata file that separates genomes into groups to compare, and will run some basic stats to determine the enrichment of modules or pfam clans between and within the groups. 

See the [enrichment help page](https://github.com/geronimp/enrichm/wiki) for more


# pathway
Pathway reads in a KO matrix and generates a Cytoscape-readable metabolic network and metadata file. Only reactions that are possible given the KOs present in the input matrix are shown, and the modules and reactions that are included in the output can be customized.

See the [pathway help page](https://github.com/geronimp/enrichm/wiki) for more


# explore
Explore is similar to pathway, but rather than generating a specified pathway it will start from a given query compound ID, and explore the possible reactions that use that compound given the enzymes present in the input KO matrix.

See the [explore help page](https://github.com/geronimp/enrichm/wiki) for more
