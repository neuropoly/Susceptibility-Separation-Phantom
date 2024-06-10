unring - tool for removal of the Gibbs ringing artefact

based on the algorithm in the publication 
Kellner, E, Dhital B., Kiselev VG and Reisert, M. 
Gibbs‐ringing artifact removal based on local subvoxel‐shifts. 
Magnetic resonance in medicine, 76(5), 1574-1581.


Content
-------

matlab - contains MATLAB/mex implementation of the algorithm, you need the fftw3.h header to compile the mex file via
	 mex ringRm.cpp -lfftw3
     this usually links to MATLAB's libfftw3 library
     In case of a compile error, the switch -compatibleArrayDims might help
         
         
fsl - contains an implementation within the FSL framework, you need fsl sources to compile, but 
      a binary is included compiled on a 64bit linux system (kubuntu)


javascript - a javascript implementation. Runs in nodejs or the webbrowser
