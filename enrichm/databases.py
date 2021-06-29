#!/usr/bin/env python3
# Imports
import os
import logging
import pickle
# Local
from enrichm.data import Data

###############################################################################

class Databases:

    def __init__(self):
        if os.path.isfile(os.path.join(Data.DATABASE_DIR, 'VERSION')):

            with open(os.path.join(Data.DATABASE_DIR, 'VERSION')) as out_io:
                self.DB_VERSION = out_io.readline().strip().replace('.tar.gz', '')

            self.CUR_DATABASE_DIR = os.path.join(Data.DATABASE_DIR, self.DB_VERSION)

            with open(os.path.join(self.CUR_DATABASE_DIR, 'VERSION')) as out_io:
                self.PICKLE_VERSION = out_io.readline().strip()

            self.IDS_DIR = os.path.join(self.CUR_DATABASE_DIR, 'ids')
            self.REF_DIR = os.path.join(self.CUR_DATABASE_DIR, 'databases')
            self.KO_HMM_CUTOFFS = os.path.join(self.CUR_DATABASE_DIR, 'ko_cutoffs.tsv')

            self.PICKLE = 'pickle'
            self.HMM_SUFFIX = '.hmm'
            self.DMND_SUFFIX = '.dmnd'
            self.KO_DB_NAME = 'uniref100.KO'
            self.EC_DB_NAME = 'uniref100.EC'
            self.PFAM_DB_NAME = 'pfam'
            self.KO_HMM_DB_NAME = 'ko'
            self.TIGRFAM_DB_NAME = 'tigrfam'
            self.CAZY_DB_NAME = 'cazy'

            self.M2DEF = os.path.join(self.CUR_DATABASE_DIR, 'module_to_definition')
            self.M = os.path.join(self.CUR_DATABASE_DIR, 'module_descriptions')
            self.COMPOUND_DESC = os.path.join(self.CUR_DATABASE_DIR, 'br08001')
            self.R2K = os.path.join(self.CUR_DATABASE_DIR, 'reaction_to_orthology')
            self.R2C = os.path.join(self.CUR_DATABASE_DIR, 'reaction_to_compound')
            self.R2M = os.path.join(self.CUR_DATABASE_DIR, 'reaction_to_module')
            self.M2R = os.path.join(self.CUR_DATABASE_DIR, 'module_to_reaction')
            self.M2C = os.path.join(self.CUR_DATABASE_DIR, 'module_to_cpd')
            self.R2P = os.path.join(self.CUR_DATABASE_DIR, 'reaction_to_pathway')
            self.P2R = os.path.join(self.CUR_DATABASE_DIR, 'pathway_to_reaction')
            self.C2R = os.path.join(self.CUR_DATABASE_DIR, 'compound_to_reaction')
            self.C = os.path.join(self.CUR_DATABASE_DIR, 'compound_descriptions')
            self.R = os.path.join(self.CUR_DATABASE_DIR, 'reaction_descriptions')
            self.P = os.path.join(self.CUR_DATABASE_DIR, 'pathway_descriptions')
            self.K = os.path.join(self.CUR_DATABASE_DIR, 'ko_descriptions')

            self.PFAM2CLAN = os.path.join(self.CUR_DATABASE_DIR, 'pfam_to_clan')
            self.PFAM2NAME = os.path.join(self.CUR_DATABASE_DIR, 'pfam_to_name')
            self.PFAM2DESCRIPTION = os.path.join(self.CUR_DATABASE_DIR, 'pfam_to_description')
            self.EC2DESCRIPTION = os.path.join(self.CUR_DATABASE_DIR, 'ec_to_description')
            self.TIGRFAM2DESCRIPTION = os.path.join(self.CUR_DATABASE_DIR, 'tigrfam_descriptions')
        else:
            raise Exception(f"\nNo database version file found. Have you: \n\
- Installed the EnrichM database using the 'enrichm data' command?\n\
- Specified the location of the EnrichM database by exporting a \
bash variable called ENRICHM_DB? (Currently I'm looking here: {Data.DATABASE_DIR})")
        self.signature_modules = set(['M00611', 'M00612', 'M00613', 'M00614',
                                      'M00617', 'M00618', 'M00615', 'M00616',
                                      'M00363', 'M00542', 'M00574', 'M00575',
                                      'M00564', 'M00660', 'M00664', 'M00625',
                                      'M00627', 'M00745', 'M00651', 'M00652',
                                      'M00704', 'M00725', 'M00726', 'M00730',
                                      'M00744', 'M00718', 'M00639', 'M00641',
                                      'M00642', 'M00643', 'M00769', 'M00649',
                                      'M00696', 'M00697', 'M00698', 'M00700',
                                      'M00702', 'M00714', 'M00705', 'M00746'])

        self.KO_DB = os.path.join(self.REF_DIR, self.KO_DB_NAME + self.DMND_SUFFIX)
        self.EC_DB = os.path.join(self.REF_DIR, self.EC_DB_NAME + self.DMND_SUFFIX)

        self.PFAM_DB = os.path.join(self.REF_DIR, self.PFAM_DB_NAME + self.HMM_SUFFIX)
        self.KO_HMM_DB = os.path.join(self.REF_DIR, self.KO_HMM_DB_NAME + self.HMM_SUFFIX)
        self.TIGRFAM_DB = os.path.join(self.REF_DIR, self.TIGRFAM_DB_NAME + self.HMM_SUFFIX)
        self.CAZY_DB = os.path.join(self.REF_DIR, self.CAZY_DB_NAME + self.HMM_SUFFIX)
        self.PFAM_CLAN_DB = os.path.join(self.IDS_DIR, 'PFAM_CLANS.txt')

    def m2def(self):
        logging.debug("Loading module descriptions")
        return self.load_pickle(self.M2DEF)

    def m(self):
        logging.debug("Loading reaction to pathway information")
        return self.load_pickle(self.M)

    def r2p(self):
        logging.debug("Loading pathway to reaction information")
        return self.load_pickle(self.R2P)

    def p2r(self):
        logging.debug("Loading reaction to orthology information")
        return self.load_pickle(self.P2R)

    def r2k(self):
        logging.debug("Loading reaction to module information")
        return self.load_pickle(self.R2K)

    def r2m(self):
        logging.debug("Loading module to reaction information")
        return self.load_pickle(self.R2M)

    def m2r(self):
        logging.debug("Loading module to compound information")
        return self.load_pickle(self.M2R)

    def r2c(self):
        logging.debug("Loading compound to reaction information")
        return self.load_pickle(self.R2C)

    def c2r(self):
        logging.debug("Loading compound descriptions")
        return self.load_pickle(self.C2R)

    def c(self):
        logging.debug("Loading pathway descriptions")
        return self.load_pickle(self.C)

    def p(self):
        logging.debug("Loading reaction descriptions")
        return self.load_pickle(self.P)

    def r(self):
        logging.debug("Loading ko descriptions")
        return self.load_pickle(self.R)

    def k(self):
        logging.debug("Loading compound classifications")
        return self.load_pickle(self.K)

    def compound_desc_dict(self):
        logging.debug("Loading pfam to clan information")
        return self.load_pickle(self.COMPOUND_DESC)

    def pfam2clan(self):
        logging.debug("Loading clan descriptions")
        return self.load_pickle(self.PFAM2CLAN)

    def pfam2description(self):
        logging.debug("Loading ec descriptions")
        return self.load_pickle(self.PFAM2DESCRIPTION)

    def ec2description(self):
        logging.debug("Loading pfam hierarchy")
        return self.load_pickle(self.EC2DESCRIPTION)

    def tigrfamdescription(self):
        logging.debug("Loading reference db paths")
        return self.load_pickle(self.TIGRFAM2DESCRIPTION)

    def k2r(self):
        k2r = dict()
        for reaction, kos in self.r2k().items():
            for ko in kos:
                if ko not in k2r:
                    k2r[ko] = list()
                k2r[ko].append(reaction)
        return k2r

    def c2m(self):
        c2m = dict()

        for module, compounds in self.m2c().items():
            substrates = compounds[0]
            for substrate in substrates:
                if substrate in c2m:
                    c2m[substrate].append(module)
                else:
                    c2m[substrate] = [module]
        return c2m

    def load_pickle(self, file):

        with open('.'.join([file, self.PICKLE_VERSION, self.PICKLE]), 'rb') as file_io:
            loaded_pickle = pickle.load(file_io)

        return loaded_pickle

    def parse_ko_cutoffs(self):
        cut_ko = dict()
        out_io = open(self.KO_HMM_CUTOFFS)
        _ = out_io.readline()

        for line in out_io:
            sline = line.strip().split('\t')

            if sline[1] == '-':
                cut_ko[sline[0]] = [0.0, "NA"]
            else:
                cut_ko[sline[0]] = [float(sline[1]), sline[2]]

        return cut_ko
