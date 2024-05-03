# AAV

Last Updated: 05/03/2024

## Disclaimer

**THIS REPOSITORY IS OBSOLETE AND WILL NO LONGER BE MAINTAINED**

This was a personal repository for analyzing PacBio long-read AAV sequencing data. The open source code is now being maintained and developed at [FormBio's LAAVA repo](https://github.com/formbio/laava). Please visit LAAVA to get the latest code base! Thanks!



## Pre-requisites

* Python 3.7
* R

Python libaries required:
* [pysam](https://anaconda.org/bioconda/pysam)

R packages required:
* ggplot2
* dplyr
* grid
* gridExtra

## Installation

You can directly download/clone the repo to use the scripts directly. 

```
$ git clone https://github.com/Magdoll/AAV.git
```

You can install the dependencies on your own or use one of the following conda-based options.

### Option1: Use Conda and install the pre-requisites individually

```
conda install -c bioconda pysam
conda install -c r ggplot2
conda install -c r dpylr
conda install -c r grid
conda install -c r gridExtra
```

### Option2: Use Conda with yml

Suppose you have [anaconda](https://docs.anaconda.com/anaconda/install/linux/) installed and the binary is in `$HOME/anaCogentPy37/bin`. You would add the binary to $PATH and create a new conda environment called `AAV.env`.

```
$ export PATH=$HOME/anaCogentPy37/bin:$PATH
$ conda env create -f AAV.conda_env.yml
$ source activate AAV.env
```

At this point the prompt should change to `(AAV.env) $`

## Usage

Please read the [AAV tutorial](https://github.com/Magdoll/AAV/wiki/Tutorial:-Analyzing-AAV-Data)

