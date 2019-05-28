# EnrichM
A set of comparative genomics utilities for large sets of metagenome assembled genomes (MAGs).

## Commands

```
                         _____            _      _     __  __ 
                        | ____|_ __  _ __(_) ___| |__ |  \/  |
                        |  _| | '_ \| '__| |/ __| '_ \| |\/| |
                        | |___| | | | |  | | (__| | | | |  | |
                        |_____|_| |_|_|  |_|\___|_| |_|_|  |_|
  ------------------------------------------------------------------------------------

  Annotation
    annotate        -> Genome annotation pipeline.

  Enrichment analysis
    classify        -> Determine what pathways a genome encodes.
    enrichment      -> Calculate enrichment of functional genes between groups.

  Network analysis
    pathway         -> Generate a metabolic network from specific KEGG module or
                       compounds.
    explore         -> Take steps into metabolism using a KEGG compound ID as a 
                       starting point. Useful to see which pathways use a compound
                       of interest.

  Machine learning
    generate        -> Generate a random forest model.
    predict         -> Run random forest model on new data.
```

## Annotate
The _annotate_ function of EnrichM is used to provide your genome with one or more of KO, PFAM, CAZy, EC, and TIGRfam annotations. KO and EC numbers are annotated by blasting against UniProt using DIAMOND, and PFAM and TIGRFAM are searched against their respective databases using hmmsearch. CAZy annotations are acquired by using the HMM database [dbCAN](http://csbl.bmb.uga.edu/dbCAN). Default cutoffs are used (see enrichm annotate -h), but these can be specified by the user. The output will include a matrix for each annotation type (IDs are rows, genomes/samples are columns) that is filled with count information. You will also be provided with a GFF file for each genome, resembling that which would be provided by prokka, but rather than string descriptions, each gene’s annotations will be provided as annotation ids (KO, PFAM/TIGRfam name). Genomes or assemblies can either be specified as:

* `--genome_files` A space separated list of paths to genomes.
* `--genome_directory` A directory of genomes. By default, EnrichM will search the directory for files ending in .fna although this can be changed with the `--suffix` flag.
* `--protein_files` A space separated of paths to amino acid sequence (from genomes) files.
* `--protein_directory` A directory of amino acid sequence files. By default, EnrichM will search the directory for files ending in .faa although this can be changed with the `--suffix` flag.

An example run to annotate a directory full of genomes (bin) using CAZy, KO, TIGRfam, EC, and Pfam annotations might be:

```
$ enrichm annotate \
	# Directiory containing genomes
	--genome_directory bin \
	# Number of genomes to run in parallel while annotating
	--threads 30 \
	# Annotate with ko ids
	--ko \
	# Annotate with pfam ids
	--pfam \
	# Annotate with tigrfam ids
	--tigrfam \
	# Annotate with cazy ids
	--cazy
    # Annotate with homologous clusters
    --cluster
    # Annotate with putative orthologs
    --orthologs
```

The output is a directory, containing the annotations and gff files for each of the genomes. By default all of the processing files of the annotation are retained, and can be found in clearly labelled sub-directories. The outputs are:

* `annotate.log` A log recording the command ran and information about the run
* `ko_frequency_table.tsv` KO table, with KO ids as rows and columns as genomes
* `pfam_frequency_table.tsv` PFAM table, with PFAM ids as rows and columns as genomes
* `tigrfam_frequency_table.tsv` TIGRfam table, with TIGRfam ids as rows and columns as genomes
* `cazy_frequency_table.tsv` CAZy table, with CAZy ids as rows and columns as genomes
* `annotations_ko` Directory containing DIAMOND BLASTp outputs from KO annotation
* `annotations_pfam` Directory containing HMMsearch outputs from PFAM annotation
* `annotations_tigrfam` Directory containing HMMsearch outputs from TIGRfam annotation
* `annotations_cazy` Directory containing HMMsearch outputs from CAZy annotation
* `annotations_gff` Directory containing gff files for each genome
* `genome_proteins` Directory .faa protein files for each genome

## Classify
The classify function is used to interpret KO annotations as KEGG modules. This is done using KEGG module definitions, which are defined on the [KEGG module page](https://www.genome.jp/kegg/module.html). This format can be used to define different series of reactions that are used to carry out a reaction. The definition can include multiple optional pathways, and EnrichM will simply check if a given genome can encode any of these pathways. Using KEGG module definition format, you can also specify your own custom KEGG modules to EnrichM in order to quickly check a large set of genomes are capable of carrying out a type of metabolism you’re interested in. For example there is no module definition for xylose degradation, in [Woodcroft et al. 2018](https://www.nature.com/articles/s41586-018-0338-1) we had to specify our own:

```
canonical_xylose_deg	K01805 K00854
fungal_xylose_deg	(K17743,K00011) (K05351,K00008) K00854
```

All you need to provide to enrichM classify is a KO matrix, which is produced by EnrichM annotate. If you have any, you can specify your custom modules within a file to the `--custom_modules` flag. The custom modules file must be a tab separated file, with the first column being the module descriptions, and the second column being the definition. You can filter the annotations EnrichM by completeness using the `--cutoff` flag, where 0 means no cutoff and 1 means only complete modules will be returned. An example run including some custom modules would look something like this:

```
$ enrichm classify --custom_modules custon_modules.tsv \
	# The input ko matrix for classification
	--genome_and_annotation_matrix ko_matrix.tsv \
	# Show only 100% complete modules
	--cutoff 1 \
	# Output file
	--output output_file  
```

## enrichment

The enrichment function is used to compare user-defined sets of genomes with each-other. an annotate output from EnrichM annotate is required to run `enrichm enrichment`. Two tests are run to compare groups, an enrichment test and an over representation test. The enrichment test looks at the presence/absence of annotations within one group vs the other and tests if annotations of a particular type are \textbf{enriched} within one of the sets. The overrepresentation test takes into account the number of each annotation present, because although an annotation may be present in two groups, one group may encode significantly more of that annotation (e.g. genome group 1 encodes far more P450s than group 2). The tests that are run to calculate enrichment are decided using the following decision tree:

```
group_1 # Number of genomes in group 1
group_2	# Number of genomes in group 2

if(group_1 == 1 or group_2 == 1):
	enrichment_test		= Presence/absence
	overrepresentation_test	= Z score test
else:
	enrichment_test		= Fisher's exact test
	overrepresentation_test	= Mann-Whitney u test
```

Muti-test correction is applied after the enrichment tests have been run. By default this is Benjamini-Hochberg FDR, but different methods can be selected from the different options available in the [statsmodels module](https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html) If analysing a KO matrix, the KOs found to be differentially abundant between the user-defined groups are put into pathways. An example run may look like the following:

```
enrichm enrichment \
	# The annotate output
	--annotate_output annotate_output \ 
	# The metadata file, defining which 
	# groups the genomes are in
	--metadata cyano_metadata.tsv \ 
	# GTDB genomes to compare with
	--batchfile cyano_gtdb.tsv \ 
	# Compare KO annotations
	--ko \ 
	# Define the output file
	--output enrichment_output
```

If you wish to compare with high quality genomes from the GTDB.

## Pathway
Create curated metabolic networks from your KO annotations. The network generating scripts can be used on genomes, assemblies, but can also be used to integrate metagenomic, transcriptomic and metabolomic data.

## Traverse
Explore which compounds in a metabolic network are "visited" the most. 

## Generate
Generate a random forest model of the genome form input annotations

## Predict
Run the random forest classifier produced but Generate on new data.