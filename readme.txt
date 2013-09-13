README: Code for simulating a tilting BEC with vortices.

CORE PROGRAM
The code is executed through the makefile.
Ie. go to the terminal and type "make", then run the program with 
"./BEC_evolve".

The main algorithm is contained within "BEC_evolve.c" which is linked
to its dependent files through the header file "BEC_3D.h". This header
file also contains the parameters for setting up the program such as
grid size, tilt angle, time etc.

The dependent files include "initial_cond.c" which sets up the grid,
potentials and the initial condition.

"matrix_functions.c" includes code for manipulating matrices,
integrating the wavefunction for printing and normalisation.

"operators.c" is dedicated to the linear (diffusion) and non-linear
(including the potential) operators. It also contains the function
which imprints a vortex.
I am assuming that it is here you would implement any kind of GPU
processing through the FFT algorithm.

These are brief descriptions, but I have tried to add more detail in 
each file for clarity.

"Code_Strucure.png" is an image used in the report for this project
which outlines the program structure. "Code_timeline.png" also may
help you to follow the codes chronological progression.

The algorithm and details are based on those used in Caradoc-Davies PhD thesis, which I have also included. Details on the algorithm can be found in Appendix B. I utilised the "Efficient implementation in Python" algorithm.

Appendix C details the diffusion operator and relevent coefficients for computing this term.


OUTPUT
Data is output in 3 ways into 3 files.

"outputXY.txt": The modulus of the wavefunction is computed and then integrated through the Z-axis to obtain the top down view we would expect from experimental observation. The data is printed to text file along with key parameters relevent for plotting the data.

"outputYZ.txt": The modulus of the wavefunction is computed and then integrated through the X-axis to obtain a side on view of the condensate, this is primarily for analysing how the tilt speed and ramping function affect BEC oscillations. The data is printed to text file along with key parameters relevent for plotting the data.

"phaseXY.txt": The phase of the complex valued wavefunction is computed at a set level XY plane. This is primarily for observing the vortex behaviour. The data is printed to text file along with key parameters relevent for plotting.


PLOTTING/ANIMATIONS
I use octave for plotting each frame of the BEC evolution. These print to a folder in the main directory. I print the plots as .gif files primarily because Octave seems to spit these out much quicker. Also .png/.bmp files suffered from some rendering defects which I could not avoid.

The .gif's are then merged into a .gif animation using the program "gifsicle".

This .gif is then converted to an .avi movie using "ffmpeg".

It is a convoluted process, but I tried a number of ways and this was the quickest! Let me know if you find a better procedure.

"printbecXY.m" takes the data from "outputXY.txt" and plots/animates the top down view of the condensate.

"printbecYZ.m" takes the data from "outputYZ.txt" and plots/animates the side view of the condensate.

"printphase.m" takes the data from "phaseXY.txt" and plots/animates the phase map of the condensate.
