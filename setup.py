import subprocess
import sys
from setuptools import setup, find_packages

def install_dependencies():
   try:
       # Check if the required environment is installed
       subprocess.check_call(['conda', 'env', 'create', '-f', 'env_thale_BSA.yml'])

       # Environment already exists
       print("Required environment is already installed.")
   except subprocess.CalledProcessError:
       # Install the environment
       print("Installing required environment...")
       subprocess.check_call(['conda', 'env', 'update', '-f', 'env_thale_BSA.yml'])
       print("Required environment installed successfully.")

install_dependencies()


setup(
    name="thale_bsa",
    version="1.0.0",
    author="Tim Shull",
    author_email="less.tact_0m@icloud.com",
    description="A python library for bulk segregant analysis of EMS-derived point mutations in Arabidopsis thaliana",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TeaShull/PyAtBSA",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bioinformatics",
    ],
    keywords="Arabidopsis BSA Genomics",
    entry_points={
        "console_scripts": [
            "thale_bsa=src.__main__:main",
        ],
    },
)

