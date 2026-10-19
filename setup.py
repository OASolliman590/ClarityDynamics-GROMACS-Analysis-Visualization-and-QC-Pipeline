#!/usr/bin/env python3
"""
Setup script for GROMACS Analysis Toolkit
"""

from setuptools import setup, find_packages
import os

# Read README for long description
def read_file(filename):
    with open(os.path.join(os.path.dirname(__file__), filename), encoding='utf-8') as f:
        return f.read()

setup(
    name='gromacs-analysis-toolkit',
    version='1.0.0',
    description='Unified toolkit for GROMACS molecular dynamics analysis and visualization',
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    author='Your Team',
    author_email='your.email@example.com',
    url='https://github.com/yourusername/gromacs-analysis-toolkit',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy>=1.20.0',
        'pandas>=1.3.0',
        'plotly>=5.0.0',
        'kaleido>=0.2.1',
        'scipy>=1.7.0',
        'pyyaml>=6.0',
        'matplotlib>=3.4.0',
    ],
    extras_require={
        'trajectory': [
            'MDAnalysis>=2.0.0',
        ],
        'prolif': [
            'MDAnalysis>=2.0.0',
            'prolif>=2.0.0',
            'rdkit>=2022.9.1',
            'ipython',  # ProLIF ligand network rendering
        ],
        'clustering': [
            'scikit-learn>=1.0.0',
        ],
        'analysis': [
            'seaborn>=0.11.0',
        ],
        'all': [
            'MDAnalysis>=2.0.0',
            'prolif>=2.0.0',
            'rdkit>=2022.9.1',
            'ipython',
            'scikit-learn>=1.0.0',
            'seaborn>=0.11.0',
        ],
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=21.0',
            'flake8>=3.9',
            'mypy>=0.910',
        ],
    },
    entry_points={
        'console_scripts': [
            'gromacs-md=gromacs_analysis.cli.md_cli:main',
            'gromacs-mmpbsa=gromacs_analysis.cli.mmpbsa_cli:main',
            'gromacs-pca=gromacs_analysis.cli.pca_cli:main',
            'gromacs-pipeline=gromacs_analysis.cli.pipeline_cli:main',
            'gromacs-orchestrate=gromacs_analysis.cli.orchestrate_cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    python_requires='>=3.8',
    keywords='gromacs molecular-dynamics visualization mmpbsa pca analysis',
)
