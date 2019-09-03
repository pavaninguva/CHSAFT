# CH-SAFT

## Introduction

This project is designed to test various equations of states including SAFT in understanding the thermodynamics of polymer blends and how immiscible blends undergo demixing using a modified Cahn-Hilliard framework. 

## Installation

### Cahn-Hilliard 

There are two ways to set up the environment with relevant dependencies for the CH side of the project. The first is using Anaconda environments, and the second is to use a Singularity container. The conda env is better suited for local development as Singularity cannot be run on a Windows WSL or Mac natively. The Singularity container is better suited for HPC use. 

Essentially, the two main dependencies are Fenics >= 2018.1.0 and pyyaml. 

#### Conda

To set up the conda env: 

```bash
conda create -n fenicsproject -c conda-forge fenics
```

You also need to install the pyyaml within the environment

```bash
conda activate fenicsproject
pip install pyyaml 
```

Another alternative is to create the environment with the available `environment.yml` file in the repo: 

```bash
conda env create -f ./fenics2019/environment.yml
```

#### Singularity

To build a Singularity container `foo.sif` with the relevant dependencies: 

```bash
singularity build foo.sif ./fenics2019/test.def
```

If you are on a cluster or do not have sudo access, an alternative might be to pull a container with the relevant dependencies: 

```bash
singularity pull library://trlandet/default/ocellaris:2019.0.2
```



## Running simulations

### Cahn-Hilliard

The main code is captured in `binary.py`. To run the script by itself: 

Within the conda env:

```bash
conda activative fenicsproject
python binary.py
```

Using singularity, assuming the container is in the current directory: 

```bash
singularity exec ./ocllaris:2019.0.2.sif python3 binary.py
```

To change various parameters in the simulation, you need to modify the `params.yml` file:

```yaml
# mole fraction of A
A_RAW: 0.5

# Thermo model
GIBBS: FH

# Materials
MATERIAL_CHOICE: PS_PMMA
...
```

The various thermodynamic models are captured within the binary.py script and can be selected based on selecting an appropriate model ```GIBBS``` in the `params.yml` file. For instance, `FH` corresponds to the Flory-Huggins model. 