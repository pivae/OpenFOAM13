#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=64 #number of processes needed
#SBATCH --exclusive #no other jobs allowed on the node while this job is running
#SBATCH --time=100:00:00
#SBATCH -J D_largeVessel #job name
#SBATCH --mail-type=ALL
#SBATCH --account=ear2127
#SBATCH --partition=Rome
#SBATCH --mail-user=enrico.pivato@vki.ac.be
#SBATCH --mem=800G
#SBATCH -o /data/nobackup/EAR2127/pivato/ARTIV/pumpsToInletCoreLargerVesselURANS_highSwirl/output.%J
#SBATCH -e /data/nobackup/EAR2127/pivato/ARTIV/pumpsToInletCoreLargerVesselURANS_highSwirl/error.%J
 
cd /data/nobackup/EAR2127/pivato/ARTIV/pumpsToInletCoreLargerVesselURANS_highSwirl

singularity exec --cleanenv -B /data/nobackup/EAR2127/pivato/ARTIV/pumpsToInletCoreLargerVesselURANS_highSwirl /data/nobackup/EAR2127/pivato/TEST/OpenFOAMv12_PCTD1.sif /bin/bash -c '
export PATH=/software/alternate/fk/openmpi/4.0.2/bin:$PATH &&
export LD_LIBRARY_PATH=/software/alternate/fk/openmpi/4.0.2/lib:$LD_LIBRARY_PATH &&
source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc &&
export FOAM_USER_LIBBIN=/opt/openfoam-12-vki/platforms/linux64GccDPInt32Opt/lib &&
export LD_LIBRARY_PATH=/opt/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/lib:$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH &&
. $WM_PROJECT_DIR/bin/tools/RunFunctions && runApplication reconstructPar -latestTime'

 

