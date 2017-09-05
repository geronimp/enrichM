# Downloads the necessary databases required to run enrichM:
# 		- PFAM
# 		- TIGRFAM
# 		- KO annotated portion of UniProt100


# Constants
DATADIR=`echo data/databases`

# Setup directory
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Creating database directory...
mkdir $DATADIR

# Download TIGRFAMs
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Downloading TIGRFAMs...
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz -O $DATADIR\/tigrfam.hmm.tar.gz
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Decompressing TIGRFAM HMMs...
tar -xvzf $DATADIR\/tigrfam.hmm.tar.gz > /dev/null
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Concatenating TIGRFAM HMMs...
cat TIGR*HMM > data/databases/tigrfam.hmm
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Cleaning up...
rm TIGR*HMM
rm data/databases/tigrfam.hmm.tar.gz

# Download PFAM
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Downloading PFAM...
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz -O $DATADIR\/pfam.hmm.gz
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Decompressing PFAM HMMs...
gunzip $DATADIR\/pfam.hmm.gz
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Cleaning up...

#Download Uniref100
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Downloading UniRef100 sequences...
wget https://data.ace.uq.edu.au/public/enrichm/data/uniref100.dmnd -O $DATADIR\/uniref100.dmnd
echo [`date '+%Y-%m-%d %H:%M:%S'` - ENRICHM SETUP] Done! Cleaning up...
