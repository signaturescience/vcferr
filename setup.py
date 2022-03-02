import pathlib
from setuptools import setup

HERE=pathlib.Path(__file__).resolve().parent
README=(HERE/"README.md").read_text()

setup(
	name='vcferr',
	version='1.0.0',
	description='Probabilistic VCF genotype error simulation',
	long_description=README,
	long_description_content_type="text/markdown",
	url='https://github.com/signaturescience/vcferr',
	author='VP Nagraj',
	install_requires=['pysam'],
	author_email='',
	entry_points='''
        	[console_scripts]
        	vcferr=vcferr.__main__:vcferr    
		''',
	packages=setuptools.find_packages(),
	zip_safe=False)
