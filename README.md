Instructions for Running PyICM for an Idealized Estuary
J. Blake Clark Sep 2021

pyICM is a simplified version of the estuarine and coastal mechanistic model developed by the Army Corps of Engineers for Chesapeake Bay, called CE-QUAL-ICM (Cerco and Cole, 1993). It can be used to better understand estuarine water quality and biogeochemistry and to test new formulations. It has been modified to calculating the fully hyperspectral underwater light field (Clark et al. 2020) using spectral data from Chesapeake Bay collected in the Rhode River, MD (Rose et al. 2018). Clark et al. 2020 contains a full model description in the Appendices that are used in pyICM
1)	Cerco, C. F., & Cole, T. (1993). Three-dimensional eutrophication model of Chesapeake Bay. Journal of Environmental Engineering, 119(6), 1006–1025.

2)	Clark, J. B., Long, W., & Hood, R. R. (2020). A comprehensive estuarine dissolved organic carbon budget using an enhanced biogeochemical model. Journal of Geophysical Research. https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JG005442

3)	Rose, K. C., Neale, P. J., Tzortziou, M., Gallegos, C. L., & Jordan, T. E. (2018). Patterns of spectral, spatial, and long-term variability in light attenuation in an optically complex sub-estuary. Limnology and Oceanography. https://aslopubs.onlinelibrary.wiley.com/doi/abs/10.1002/lno.11005

This software has been tested on macOS Catalina version 10.15.7.
Before we start
Open a terminal by searching in the finder for “Terminal” and opening the window
The terminal is the primary way to install and run the programs for pyICM. In a future release it will be packaged with a simple graphical user interface but until then the terminal is necessary. These instructions are geared for Mac or Linux (i.e. Unix) based OS’s so apologies for the lack of specific information for those that use Windows. In the near future I will test and update for Windows.
Having a code editor to edit the software easily
Visual Studio Code is the best open-source code visualization and editing software which is free to download at https://code.visualstudio.com/
Spyder is the other tool that can be useful to see the variables in the workstation.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1)	Install Anaconda
a)	https://docs.anaconda.com/anaconda/install/ and follow instructions for your specific OS and version of pythoni,ii
i)	Python comes pre-installed on Mac and most Linux OS’s, but may need to check version to make sure you install the correct Anaconda package
(1)	To check version type into the terminal window python --version
ii)	Windows does not come with python pre-installed, so navigate to the link https://docs.python.org/3/using/windows.html to install the most up-to-date version of python
b)	Once installed can type conda help to see what commands are available
c)	Type conda update to make sure that everything is up to date (even if you just installed it)
d)	Type conda list to list the packages that are currently installed with the base installation
2)	Install Jupyter Notebooks in python
a)	Jupyter notebooks are super handy for visualizing code and stepping through blocks to better understand the code. There is a notebook file that comes with pyICM that can be visualized in the notebook application
b)	To install Jupyter type conda install -c conda-forge jupyterlab https://jupyter.org/install
3)	Make sure packages are installed that are needed for pyICM
a)	Most of the packages need are standard with the conda install, to see all packages open the ICM_main.py file. At the top, the packages are visible. Many of them are files that are contained within the code directory and are specific to pyICM. A brief description of the other packages are below, along with the current (as of Sep 2021) list from ICM_main.py











•	The following four come with python/conda and don’t need to be installed
o	sys system-specific parameters and functions https://docs.python.org/3/library/sys.html
o	os which is used to manipulate and change paths and directories
o	numpy used for many numerical calculations and operators https://numpy.org/
o	matplotlib used for plotting things https://matplotlib.org
o	pandas used for data storage and manipulation https://pandas.pydata.org/
•	The rest (Density - Oxygen) are part of the source code of pyICM and descriptions of each file can be found at the top in the source code for each. A brief description is provided below
o	Density calculates the density of water based on temperature and salinity
o	functions contains functions that are needed in multiple modules of pyICM
o	DefVars defines global parameters that are used across pyICM
o	Light calculates formulations for the prediction and propagation of light through the water column
o	Sediment calculates sediment sinking and transport
o	Nutrients calculates reactions of nutrients (NH4+, NO3-, and PO43-)
o	DOM calculates reactions of dissolved organic matter (dissolved organic carbon and nitrogen)
o	POM calculates reactions of particulate organic matter (particulate organic carbon and nitrogen)
o	Algae calculates reactions of phytoplankton (2 groups, cold and warm seasons)
o	Oxygen calculates reactions of dissolved oxygen
•	Additional packages that are needed (in the light module)
o	Pysolar used to calculate the position of the sun based on latitude, longitude, and date https://pysolar.readthedocs.io/en/latest/
	To install in the terminal type pip install pysolar
o	Datetime used to manipulate dates
	To install in the terminal type pip install datetime
4)	Now that all software and packages are installed, we can run the model in its initial configuration
5)	The main program is called ICM_main.py  which is where the other modules are called from, the forcing files are read in, and the model integrates in time and space before plotting some variables. I highly recommend opening the MAIN_PROGRAM jupyter file and stepping through so you can see what each block of code is doing. *MAKE SURE JUPYTER LAB IS INSTALLED*
6)	Within the main program there is a block of code where the station data is read in as the forcing and the weather from the North American Regional Reanalysis is found.  This is where initially things can be changed to simulate potential environmental changes that would come with different management strategies or climate change scenarios.
7)	The file DefVars.py is where a lot of parameters are stored that are used in the other modules. A more advanced experiment could change some of these parameters to alter how nutrients, phytoplankton, and oxygen change and interact through time. All parameters have a short explanation within this file.

	


