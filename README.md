# Parallelizing Protein Docking Algorithm using CUDA 

> # About FTDock   
>
>> The [FTDock](http://www.sbg.bio.ic.ac.uk/docking/index.html) programme is intended to be able to dock two proteins. This means starting
from the known structures of two protein subunits of a biological complex known to exist, in
unbound conformations, and ending up with a limited set of possible models for the complex.
>>
>> The ftdock algorithm is based on that of Katchalski-Katzir. It discretises the two
molecules onto orthogonal grids and performs a global scan of translational and rotational
space. In order to scan rotational space it is necessary to rediscretise one of the molecules
(for speed the smaller) for each rotation. The scoring method is primarily a surface complementarity score between the two grids, and this is shown in Figure below. To speed up the surface
complementarity calculations, which are convolutions of two grids, Fourier Transforms are
used. This means that the convolutions are replaced with multiplications in Fourier space,
and despite having to perform the forward and reverse Fourier Transforms, this decreases
the overall computation required. The surface complementarity was the only score used in
the original method. The original work on ftdock by Gabb found it a useful addition to
include an electrostatic filter, and this is again implemented in the current version. 

# History

- The original FTDOCK programme was written in 1997 by Henry Gabb. It was written in Fortran 77 with parallel capabilities for a Silicon Graphics Challenge machine. It was used for the CASP2 conference with some success.

- The Last available version of FTDock (3D-Dock suite) was written by Gidon Moont. It is written in C, with no parallel capabilities. 

- This version of FTDOCK is written in CUDA C/C++ parallelizing many functions of the programme.
</br></br>

[For scope of parallelization in this programme](https://slides.com/adityaranjanjha/code/fullscreen)

# Getting Started

Clone this repository:

    git clone https://github.com/Adiboy3112/CCDS-Research-Project.git


Change directory:

    cd CCDS-Research-Project/3D_Dock/progs/

Now to build the binaries:

    make

Now for testing add the mobile and static pdb files(Parsed) in the same directory and run the ftdock binary:

    ./ftdock -static <static_file_name>.parsed -mobile <mobile_file_name>.parsed > output



**Note** - The programme has perl scripts to convert pdb files into parsed files Refer [manual](./3D_Dock/Reference/3d-dock-manual.pdf) to run those scripts
</br></br>

Requirements:

- GNU/Linux OS or Windows WSL 
- Nvidia Cuda-Toolkit
- Nvidia GPU with compute capablity > 3.5 (This version is compiled for compute_capablity 7.5 and greater)

# Results

### Speed 

This version of ftdock completes within 5 mins for docking the bovine pancreatic trypsin inhibitor bound to kallikrein A complex where as the original code takes more than 5 hours to complete on my device  

### Accuracy

To be updated soon.......

# System Specs

This programme was tested on my pc and Param Shakti IIT KGP's Supercomputer:

## My PC Specs

CPU:  `Intel(R) Core(TM) i5-10300H`

GPU:  `Nvidia GeForce GTX 1650 Ti`

## Param Shakti Specs Used

CPU: `2* Intel Xeon SKL G-6148 per node`

GPU: `2*Nvidia V100 per node`

Partition used: `GPU-LOW`
     
