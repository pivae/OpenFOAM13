# OpenFOAM13

Repository containing OpenFOAM 13 material developed during my work at VKI.

## Contents

- `finiteVolume/`
  
  Includes an implementation of the **PCTD correlation**, which is not available in the standard OpenFOAM 13 distribution.

- `OpenFOAM13_loading.sh`
  
  Simple utility script to source within a SLURM job script.  
  It automatically loads the correct OpenFOAM environment depending on the cluster partition being used.

## Notes

The repository is mainly intended for:
- custom OpenFOAM developments,
- cluster environment setup,
- and VKI-related CFD workflows.
