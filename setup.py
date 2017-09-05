from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('enrichm/version.py').read()) # loads __version__

setup(name='enrichm',
      version=__version__,
      author='Joel Boyd, Ben Woodcroft, Alex Baker',
      author_email='joel.boyd@uqconnect.edu.au',
      description='enrichm is a toolbox for comparing the functional composition of population genomes',
      long_description=readme,
      license='GPL3+',
      keywords=["MAGs", "Population genomes", "metagenomics", "Annotation", "Comparison"],
      packages=find_packages(exclude='docs'),
      install_requires=('python-dateutil >=2.5.1',
                        'statsmodels >=0.8.0rc1',
                        'numpy >=1.9.1',
                        'pandas >=0.17.1',
                        'scipy >=0.17.0',
                        'biopython >=1.66'),
      setup_requires=['nose>=1.0'],
      test_suite='nose.collector',
      url='https://github.com/geronimp/enrichM',
      download_url='https://github.com/geronimp/enrichm/archive/%s.tar.gz' % (__version__),
      scripts=['bin/enrichm'],
      data_files=[
          ('pickles', ['data/br08001.07-08-2017.pickle',
                       'data/clan_to_name.07-08-2017.pickle',
                       'data/clan_to_pfam.07-08-2017.pickle',
                       'data/compound_descriptions.07-08-2017.pickle',
                       'data/compound_to_reaction.07-08-2017.pickle',
                       'data/module_descriptions.07-08-2017.pickle',
                       'data/module_to_definition.07-08-2017.pickle',
                       'data/module_to_reaction.07-08-2017.pickle',
                       'data/pathway_descriptions.07-08-2017.pickle',
                       'data/pathway_to_reaction.07-08-2017.pickle',
                       'data/pfam_to_clan.07-08-2017.pickle',
                       'data/pfam_to_description.07-08-2017.pickle',
                       'data/pfam_to_name.07-08-2017.pickle',
                       'data/reaction_descriptions.07-08-2017.pickle',
                       'data/reaction_to_compound.07-08-2017.pickle',
                       'data/reaction_to_module.07-08-2017.pickle',
                       'data/reaction_to_orthology.07-08-2017.pickle',
                       'data/reaction_to_pathway.07-08-2017.pickle'],

           'ids',     ['data/ids/KO_IDS.txt',
                       'data/ids/PFAM_CLANS.txt',
                       'data/ids/PFAM_IDS.txt',
                       'data/ids/TIGRFAM_IDS.txt']),
      ],
)
