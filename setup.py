from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('enrichm/version.py').read())  # loads __version__

setup(
    name='enrichm',
    version=__version__,
    author='Joel Boyd, Ben Woodcroft, Alexander Baker',
    author_email='joel.boyd@uqconnect.edu.au',
    description='enrichm is a toolbox for comparing the functional composition of population genomes',
    long_description=readme,
    long_description_content_type='text/markdown',
    license='GPL3+',
    keywords=['MAGs', 'Population genomes', 'metagenomics', 'Annotation', 'Comparison'],
    url='https://github.com/geronimp/enrichM',
    packages=find_packages(exclude=['docs', 'test']),
    include_package_data=True,
    scripts=['bin/enrichm'],
    python_requires='>=3.9',
    install_requires=[
        'python-dateutil>=2.8.0',
        'statsmodels>=0.14.0',
        'numpy>=1.24',
        'pandas>=2.0',
        'scipy>=1.10',
        'biopython>=1.80',
        'fuzzywuzzy>=0.18.0',
        'six>=1.16.0',
        'scikit-learn>=1.3',
        'decorator>=5.0',
        'pyarrow>=14.0',
        'certifi>=2024.0',
        'polars>=1.0',
        'pyrodigal>=3.0',
    ],
)
