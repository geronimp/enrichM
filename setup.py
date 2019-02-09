from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('enrichm/version.py').read()) # loads __version__

setup(name='enrichm',
      version=__version__,
      author='Joel Boyd, Ben Woodcroft, Alexander Baker',
      author_email='joel.boyd@uqconnect.edu.au',
      description='enrichm is a toolbox for comparing the functional composition of population genomes',
      long_description=readme,
      license='GPL3+',
      keywords=["MAGs", "Population genomes", "metagenomics", "Annotation", "Comparison"],
      packages=find_packages(exclude='docs'),
      install_requires=('python-dateutil>=2.5.1',
                        'statsmodels>=0.8.0rc1',
                        'numpy>=1.9.1',
                        'pandas>=0.17.1',
                        'scipy>=0.17.0',
                        'biopython>=1.66',
                        'six>=1.10.0',
                        'scikit-learn>=0.0'),
      setup_requires=['nose>=1.0'],
      test_suite='nose.collector',
      url='https://github.com/geronimp/enrichM',
      download_url='https://github.com/geronimp/enrichm/archive/%s.tar.gz' % (__version__),
      scripts=['bin/enrichm']
)





