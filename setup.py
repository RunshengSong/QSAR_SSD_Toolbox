from distutils.core import setup
from setuptools import find_packages

setup(
    name='QSAR_SSD_Toolbox',
    version='0.12',
    author="Runsheng Song",
    author_email="runsheng@umail.ucsb.edu",
    packages=find_packages(),
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    include_package_data=True,
    python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
)
