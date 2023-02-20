from setuptools import setup, find_packages
from CGDesigner.__init__ import __version__, __author__, __email__, __name__

def read(path):
    with open(path, "r") as f:
        out = f.read()
    return out

long_description = read("README.md")

setup(
    name=__name__,
    version=__version__,
    license="MIT",
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,
    description="Python package for the design and analysis of conditional guide RNAs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    zip_safe=False,
    packages=find_packages(include=[__name__]),
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=open("requirements.txt").read().strip().split("\n"),
    setup_requires=open("requirements.txt").read().strip().split("\n"),
    url="https://github.com/zchen15/CGDesigner",
    keywords="cgRNA, CRISPR, conditional guide RNA, nupack, RNA analysis",
    entry_points={
        "console_scripts": [__name__+'='+__name__+'.main:main'],
    },
    classifiers=[
        "Environment :: Console",
        "Framework :: Jupyter",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Utilities",
    ],
)