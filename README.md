
******************************************************************************************
MUSIC Hydro requires an initial energy density profile in (x,y,eta_s) grids.
To run  
====================================================================
||                                                                                                                                                      ||
|| gfortran main.f90 ampt2hydro1.f90 ampt2hydro2.f90 ampteigen.f90 spline.f90 -llapack -lblas  ||
||                                                                                                                                                      ||
====================================================================
*******************!!Lapack library is needed!! sudo apt-get install *********************


In the "input_file.dat" there are 12 entries:
    1) input_method
    2) k 
    3) dx
    4) tau_ampt
    5) sigma
    6) netas
    7) detas
    8) etas_min 
    9) etas_wind (works for input_method=1 &2)
  10) sigma_etas (works for input_method=3)
  11) max_etas  (works for input_method=3)
  12) ks_0 
The details are given below:

1)    There are 3 methods to initialize :
	if (input_method == 1) ---> 2d+plateau
	if (input_method == 2) ---> 3d_method1
	if (input_method == 3) ---> 3d_method2

	for details please see document.pdf

	A. In case of 2d+plateau the smearing along eta_s (f(eta_s)) is done in "init.cpp" subroutine. Therefore, we generate initial condition at eta_s=0 only.
	     parameters corresponding to f(eta_s), i.e.  eta_fall and eta_flat are chosen from the music input file.
	B. In case of 3d_method1 distributes energy in x,y,eta_s grid.
	C. In case of 3d_method2 distributes energy in x,y,eta_s grid.
	
2)     k : determines no. of spatial grid points.
3)     dx : spatial width.

        i.e.   spatial range : -(k-1)*dx/2 to +(k-1)*dx/2
        
4)     tau_ampt : initial informations obtained from AMPT at this tau.
5)     sigma : guassian smearing parameter in x-y plane.
6)     netas :  no. of eta_s points.
7)     detas :  width alond eta_s direction.
        i.e.  eta_s range : -(netas/2)*detas to (netas/2)*detas-detas
        
8)     etas_min : 
        for input_method ==1 : all contributions from eta_s<|etas_min| are taken, thus etas_wind=2*|etas_min| .
        for input_method ==2 : etas_min= (-(netas/2)*detas)-etas_wind/2
                                              for the first point -(netas/2)*detas the contributions are taken in the window [(-(netas/2)*detas)-etas_wind/2,(-(netas/2)*detas)+etas_wind/2]
        for input_method ==3 : etas_min= (-(netas/2)*detas)            
                              
9)    etas_wind (works for input_method=1 &2 ) :     
       for input_method ==1 :  etas_wind=2*|etas_min| . 
       for input_method ==2 :  etas_wind can be chosen any number.  
         
                                                
10)  sigma_etas (works for input_method = 3 ) : guassian smearing parameter along eta_s direction.
11)  max_etas (works for input_method = 3 ) : all contributions in the window eta_s < |max_etas| are taken.
12)  ks_0 : is the normalzation factor that ensures final multiplicity.

for 3 input methods there are 3 init.cpp and music input files.

please change inputs in music input file (in my case "music_input_mode_1_myinput_preliminary") accordingly.
P.S. 1. be assure that "initialize_with_entropy" inside music input file is 0.
       2. please give the correct path of "music_input.dat" in the "music_input_mode_1_myinput_preliminary".
