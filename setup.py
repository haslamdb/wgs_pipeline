#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="wgs_pipeline",
    version="0.1.0",
    description="Whole Genome Sequencing Analysis Pipeline",
    author="Your Name",
    author_email="your.email@example.com",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "scikit-allel>=1.3.5",
        "biopython>=1.80",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
    ],
    entry_points={
        'console_scripts': [
            'wgs_pipeline=wgs_pipeline:main',
            'vcf2snpmatrix=vcf2snp_matrix:main',
        ],
    },
    python_requires='>=3.8',
)
