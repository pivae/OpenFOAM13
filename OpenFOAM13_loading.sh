#!/bin/bash

module purge
module load openmpi

# --- Select correct OpenFOAM environment based on partition ---
if [[ "$SLURM_JOB_PARTITION" == "Rome" ]]; then
    echo "Running on Rome (Zen2)"
    source /software/alternate/opensource/spack/fc30/1.1.1/zen2/spack-1.1.1/opt/spack/linux-zen2/openfoam-org-13-pnesuctqbflcz3vjckxflhab5hotmv2i/etc/bashrc

elif [[ "$SLURM_JOB_PARTITION" == "Naples" ]]; then
    echo "Running on Naples (Zen1)"
    source /software/alternate/opensource/spack/fc30/1.1.1/spack-1.1.1/opt/spack/linux-zen/openfoam-org-13-ibgio47beoqntayqqepgo3de7mj4pwb3/etc/bashrc

else
    echo "Unknown partition: $SLURM_JOB_PARTITION"
    exit 1
fi

