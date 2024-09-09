from setuptools import setup, find_packages
from src import version
__version__ = version.__version__
###


setup(
    name='cgphylo',
    version=__version__,
    packages=find_packages(),
    url='https://github.com/genomicepidemiology/cgphylo',
    license='',
    install_requires=[],
    author='Malte B. Hallgren',
    scripts=['bin/cgphylo'],
    author_email='malhal@food.dtu.dk',
    description='cgPhylo - core genome phylogenetics'
)
