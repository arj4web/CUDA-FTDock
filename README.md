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

# Work done

* Downloaded the code base from [FTDock Download](http://www.sbg.bio.ic.ac.uk/docking/download.html) page and reviewed the main code file [ftdock.c](./3D_Dock/progs/ftdock.c) along with the other reference files to understand the data structures.

* Now I converted the file to [ftdock.cu](./3D_Dock/progs/ftdock.cu) a cuda runtime c code.

* There I changed the fftw library code to cufft library code and finally removed the fftw library path from header file [structures.h](./3D_Dock/progs/structures.h) making the [ftdock.cu](./3D_Dock/progs/ftdock.cu) code error free without any syntax error and compilation error due to the fourier transformation part in main loop.


[For scope of CUDA in this programme](https://slides.com/adityaranjanjha/code/fullscreen)
