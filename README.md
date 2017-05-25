[logo](logo/logo.png)

The purpose of enrichM is threefold: First, to provide a very basic pipeline for quickly annotating population genomes. Second, to compare user-defined groups of genomes and identify the enrichment of genes or metabolic pathways within and between groups. Third, to construct metabolic networks from the annotated population genomes. enrichM is not restricted to population genomes, and any of the above functions can be used to analyze and integrate metagenomic and metatranscriptomic data in metabolic networks as well. 

enrichM is currently in a state of flux and currenty does not have unit tests so while you're welcome to use it, do so at your own risk.

enrichM currently has the following sub-commands:

'''
  Genome annotation
    annotate    -> Annotate genomes

  Enrichment analysis
    classify    -> Interpret annotations provided to a genome
    build       -> Determine what KEGG modules are encoded within a population genome.
    enrichment  -> Generate an enrichment matrix from modules produced by annotate.
    
  Network analysis
    pathway     -> Generate a metabolic network from specific KEGG module or compounds
    explore     -> Explore a metabolic network from a given compound            
'''

# annotate
Annotate is a function that allows you to annotate your population genomes with a range of databases, such as KO, PFAM, TIGRFAM, and COG. The result will be a .gff file for each genome, and a frequency matrix for each annotation type where the rows are annotation IDs and the columns are genomes. This pipeline can be applied to metagenomes

See the [annotate help page](https://github.com/geronimp/enrichm/wiki) for more

# matrix
Matrix is a ...


#classify
Classify is a ...

#build
Build is a ...

#enrichment
Enrichment is a ...

#module_ab
Module_ab is a ...

#pathway
Pathway is a ...

#explore
Explore is a ...
