#!/usr/bin/env bash

# Constants
DATADIR=`data/databases`

# Setup directory
mkdir $DATADIR

# Download PFAM
echo "Downloading PFAM"
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz -O $DATADIR\/pfam.hmm
echo "Decompressing PFAM HMMs"
gunzip $DATADIR\/pfam.hmm

# Download 
