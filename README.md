# OpenFOAM13

Repository containing OpenFOAM 13 material developed during my work at VKI.

## Contents

- `finiteVolume/`
  
  Includes an implementation of the **PCTD correlation**, which is not available in the standard OpenFOAM 13 distribution.

- `OpenFOAM13_loading.sh`
  
  Simple utility script to source within a SLURM job script.  
  It automatically loads the correct OpenFOAM environment depending on the cluster partition being used.

- `OpenFOAMv12_PCTD1.def`

  Definition to build singularity container. All internal libraries (OpenMPI, UCX, GCC, ...) are optimized to be run in VKI cluster.
  Performances similar to system installations (90% of performance with respect to system installation).

- `scriptPCTD.sh`

  Script ad hoc to properly load the libraries inside the container and being able to run **OpenFOAMv12_PCTD.sif** inside the VKI cluster. 

## Notes

The repository is mainly intended for:
- custom OpenFOAM developments,
- cluster environment setup,
- and VKI-related CFD workflows.
