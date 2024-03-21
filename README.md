# LAAVA: Long-read AAV Analysis


## Installation

You can directly download/clone the repo to use the scripts directly. 

```
$ git clone https://github.com/formbio/AAV.git
```

You can install the dependencies on your own or use one of the following options.

### Option 1: Conda

The conda (or mamba) channels and dependencies are in the configuration file `AAV.conda_env.yml`. 

First, install conda via Miniconda or Anaconda. Then, for example, suppose you have [anaconda](https://docs.anaconda.com/anaconda/install/linux/) installed and the binary is in `$HOME/anaCogentPy37/bin`. To make the installed scripts available in your environment, you would add the binary to $PATH if it isn't there already:

```
$ export PATH=$HOME/anaCogentPy37/bin:$PATH
```

Next, create a new conda environment called `AAV.env` and install these dependencies from the YAML configuration file:

```
$ conda env create -f AAV.conda_env.yml
```

Finally, once installation completes, activate the new environment:

```
$ source activate AAV.env
```

At this point the prompt should change to `(AAV.env) $` and the executable AAV scripts should be available in your PATH.


### Option 2: Docker

The `Dockerfile` in this repo installs the scripts and all their dependencies into a Docker container image that you can then use to run the scripts.

To build the container image with the name `aavqc` (you can use another name if you prefer):

```
docker build -t aavqc:latest -f Dockerfile .
```

To run the container in the current working directory:

```
docker run -v $(pwd):$(pwd) -w $(pwd) -it aavqc:latest bash
```

This opens a Bash shell with the scripts in the PATH, and the original working directory mounted in place.


### Option 3: Manual installation

The prerequisites to run these scripts include:

* Python 3.7 or later
* R 3.6 or later

Python packages:
* [biopython](https://anaconda.org/bioconda/biopython)
* [pysam](https://anaconda.org/bioconda/pysam)
* [parasail-python](https://anaconda.org/bioconda/parasail-python)

R packages:
* tidyverse
* flextable
* Rmarkdown


## Testing

The `test/` subdirectory in this repo contains a Makefile that can fetch example PacBio datasets from a public server and run the scripts on them to reanalyze them and produce example HTML and PDF reports.

Once you've completed installation (above), activate your conda environment or Docker container and change to the test directory:

```
cd test
```

To fetch the public datasets and run the complete analysis (this takes some time):

```
make slow
```

Or, to just generate the HTML and PDF reports from the existing intermediate files (this takes about 2 minutes):

```
make fast
```


## Usage

Please read the [AAV tutorial](https://github.com/Magdoll/AAV/wiki/Tutorial:-Analyzing-AAV-Data)

