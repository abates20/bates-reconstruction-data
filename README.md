# Impact of data selection on automated reconstruction of genome-scale metabolic models

This repository contains the code that was used to run the analyses in the
paper, "Impact of data selection on automated reconstruction of genome-scale
metabolic models."

## Repository structure

The code in this repository is split into three different locations:

- [**src/**](./src/) contains a local python library with methods that are
  imported in other python scripts throughout the repository. It also contains
  the local implementations of Reconstructor and CarveMe.
- [**bin/**](./bin/) contains several scripts that can be reused across multiple
  analyses. Slurm scripts are stored here.
- [**analyses/**](./analyses/) contains the driver scripts and other code that
  are specific to each analysis. There is a dedicated subfolder for each
  analysis as well as a subfolder with scripts for producing each figure.

Most input data files for the analyses are stored in [**data/**](./data/), with
a some files that are specfic to individual analyses being stored in the
corresponding analysis subfolder.

Dependencies for the project are specified in
[pyproject.toml](./pyproject.toml).

## Steps to reproduce analysis

> [!IMPORTANT]
> The analyses in this project were done on a high-performance computing (HPC)
> cluster using slurm. Consequently, you need to run the code in this repository
> in a computing envirnoment where you can submit jobs using slurm.

This section describes how to run the code in this repository to reproduce the analyses and figures in the paper.

### 1. Set up a virtual environment

> [!IMPORTANT]
> All commands in this section should be run from the root directory of the
> repository (i.e., where this README file is located).

First create a virutal environment with a compatible version of python (>=3.11).
One way to do this is using the built-in venv module in python:

```sh
python -m venv .venv
```

Then install the python dependencies in the virtual environment:

```sh
source .venv/bin/activate # activates the venv
python -m pip install -e . # installs dependencies
```

> [!TIP]
> Other tools, such as [uv](https://docs.astral.sh/uv/), can also be used to
> easily manage the virtual environment and dependencies.


### 2. Additional dependencies

Next, make sure you have the DIAMOND sequence aligner (version 2.1.6 or later)
installed/available on your system (see the [DIAMOND wiki][diamond_wiki] if you
need to install it). After installation, either add the directory containing the
DIAMOND binary to your PATH or specify the path to the DIAMOND binary in a
`.env` file (see below for info about the `.env` file).

[diamond_wiki]: https://github.com/bbuchfink/diamond/wiki/2.-Installation

Additionally, make sure you have a license installed for the Gurobi solver (you
can get an [academic license](https://www.gurobi.com/academics) for free). If
the license is saved in a non-default location, you must specify the path to the
license in the `.env` file.

Finally, specify any other environment variables where you want to use a
non-default value in the `.env` file. For example, you can provide slurm with
non-default values for account and partition by specifying them in `.env`. See
the section below about local configurations for more info.

### 3. Run analysis scripts

> [!WARNING]
> Generating all the models and data in this project requires around 90 GB of
> disk space. The 04_single_misannotations analysis is the primary culprit,
> requiring around 76 GB of disk space by itself.

Each analysis has a main `run.sh` driver script that submits all the jobs for
that analysis. You can run them individually:

```sh
bash analyses/01_reconstructor_vs_carveme/run.sh
```

Or you can run them all at once:

```sh
bash analyses/runall.sh
```

Once the jobs for an analysis have finished successfully, you can generate the
corresponding figures using the `make.sh` script in the figure scripts folder:

```sh
bash analyses/06_figure_scripts/make.sh 2a # generates figure 2a
```

If the jobs for all analyses have successfully completed, then you can create
all the figures in one go:

```sh
bash analyses/06_figure_scripts/make.sh all # generates all figures
```

## Local configurations (.env file)

If any configurations need to be set for your local environment, you can create
a file named `.env` in the root of the project and specify configurations there
using the format

```
# .env file
NAME=value
```

See [.env.example](./.env.example) for example environment variables that can be
set for this repository.
