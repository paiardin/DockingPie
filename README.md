<p align="center">
<img src="https://user-images.githubusercontent.com/68196372/228539016-ec73639f-2dae-418d-9b66-07bd9dfb2425.jpeg" height="200">
</p>
<h1 align="center">A Consensus Docking plugin for PyMOL</h1>

![immagine](https://user-images.githubusercontent.com/68196372/158363199-aaaabe39-ce47-4a3e-b794-e57e088c2cce.png)


**DockingPie** is a plugin of the popular molecular graphics system PyMOL [1] that offers a versatile and user-friendly graphical user interface for assisting molecular and consensus docking analyses. 

At the current, and first release, the implemented docking programs are **Smina**, **Autodock Vina**, **RxDock** and **ADFR** [2-5]. 

Providing an easy interface to four docking programs, DockingPie is particularly suited as a platform to carry out **consensus docking and scoring analyses**.


## NEWS 

26 November 2022 - Release of DockingPie 1.2.1

- Fixed inaccessible configuration resources from GitHub (GitHub censorship)
- Now compatible with Vina 1.2.3 


## Known Issue 

"Check for updates" functionality in the CONFIGURATION TAB is currently out of use. 

## Cite

If you use DockingPie in your work, please cite:

Serena Rosignoli and Alessandro Paiardini, DockingPie: a consensus docking plugin for PyMOL, Bioinformatics, 2022, btac452, [DOI](https://doi.org/10.1093/bioinformatics/btac452 "Rosignoli et al., Bioinformatics,2022")

<img src="https://user-images.githubusercontent.com/68196372/228542114-1d45f7b0-e5a6-4862-8790-df9118f5b88f.png" height="100">

## Requirements

**Minimal requirement**: a recent version of PyMOL installed on your computer. 

DockingPie is compatible with incentive PyMOL builds distributed by [Schrodinger](https://pymol.org/2/ "Schrodinger website") (required PyMOL version >= 2.3.4) and open source builds (required PyMOL version >= 2.3.0).

DockingPie is distributed freely to the public and it has been tested and runs on Windows, macOS and Linux versions of PyMOL.


(Some incompatibilities may arise with the usage of PyMOL version 2.5.x if ‘undo’ function is enabled, which in PyMOL 2.5.2 still shows some shortcomings. Therefore, when the plugin is opened, the ‘undo’ function is automatically disabled and it is strongly suggested to keep it disabled when using the plugin.)


## Download

DockingPie plugin ZIP file: [download from here](https://github.com/paiardin/DockingPie/archive/refs/heads/main.zip "DockingPie plugin ZIP file direct download") 


## Installation 
    
DockingPie is installed, as any other PyMOL [1] plugin, via the PyMOL plugin manager:

* First download the latest version of the plugin ZIP file [here](https://github.com/paiardin/DockingPie/archive/refs/heads/main.zip  "DockingPie plugin ZIP file direct download") 

* Launch PyMOL and use the *Plugin* → *Plugin Manager* command from the main menu of PyMOL. The plugin manager window of PyMOL will open.

* Click on *Install New Plugin* and press the *Choose File…* button. Select the **DockingPie ZIP file** which you have downloaded before. 
You will be asked to give the path of the directory in which to install the plugin files. Just select the default option if you are unsure about what to do (the location of the plugin files does not make any difference when running the plugin).

## Configuration 

DockingPie, at the current and first release, integrates four different docking programs: RxDock [2], Vina [3], Smina [4] and ADFR [5]; several chemo-informatics python modules (i.e. AutoDockTools [6], Openbabel [7], sPyRMSD [8]) and other external tools like sdsorter [9]. The CONFIGURATION tab provides an easy way for the installation of the needed tools from within the plugin in two steps, as reported next.

* Configure external tools: CONFIGURATION tab → *Configure* → *Start Download* → *Finish Download*

* Install external tools: If the needed tools are not currently installed on the user’s machine, the *Install* button is enabled and it can be used to install the external components.


## How to use DockingPie

A detailed explanation on how to use DockingPie, some videos and tutorials can be found at the following links.

User's Guide: [download from here](https://github.com/paiardin/DockingPie/releases/download/versioning/DockingPie_User_Guide.pdf "DockingPie User's guide")

DockingPie GitHub Wiki: [Home](https://github.com/paiardin/DockingPie/wiki); [User's Guide](https://github.com/paiardin/DockingPie/wiki/User's-Guide); [Tutorials](https://github.com/paiardin/DockingPie/wiki/Tutorials)

Website: [Structural Bioinformatics Group at Sapienza](http://schubert.bio.uniroma1.it/)


## Tested platforms 

| PyMOL version |          Operating system          |         PyMOL source             |      
|:-------------:|:----------------------------------:|:--------------------------------:|
|     2.5.2     | Linux (Ubuntu 20.04.3 LTS), 64-bit |          Incentive               |
|     2.5.1     | Linux (Ubuntu 21.04), 64-bit       |   Incentive (Conda package)      |
|     2.5.2     | Linux (Ubuntu 18.04.2 LTS), 64-bit |          Incentive               |
|     2.5.0     | Linux (Ubuntu 21.04), 64-bit       |   Open source (Conda package)    |
|     2.5.0     | Linux (Ubuntu 20.04.3 LTS), 64-bit |   Open source (Conda package)    |
|     2.4.1     | Windows (v.10 Home), 64-bit        |          Incentive               |
|     2.5.2     | Windows (v.10 Pro), 64-bit         |          Incentive               |
|     2.5.0     | Windows (v.10 Pro), 64-bit         |  Open source (Conda package)     |
|     2.5.2     | MacOS (High Sierra v.10.13.6)      |          Incentive               |
|     2.4.0     | MacOS (High Sierra v.10.13.6)      |  Open source (Conda package)     |
|     2.5.2     | MacOS (Monterey v.12.2.1)          |          Incentive               |
|     2.5.2     | MacOS (Big Sur v.11.6.5)           |          Incentive               |
|     2.4.0     | MacOS (Big Sur v.11.6.5)           |  Open source (Conda package)     |


## DockingPie in open-source PyMOL version

DockingPie is fully compatible with open-source PyMOL version. 
However, since the “Conda Package Manager” is not integrated in open-source PyMOL, the automatic installation of the dependencies is not fully supported from within DockingPie. Thought that RxDock, sPyRMSD and OpenBabel would be manually installed by the User, we have added a section in the [User's Guide](https://github.com/paiardin/DockingPie/releases/download/versioning/DockingPie_User_Guide.pdf "DockingPie User's guide") (section 2.3) to help handling DockingPie and its dependencies on open-source PyMOL. 


## Licence

DockingPie is distributed under **GPL-3.0 License** 


## Version 

**Current Release: DockingPie 1.0.1**

### Version legend
Major changes **(1.x.x)**: a complete new version of the plugin is provided. Re-installation is needed.

Minor changes **(x.0.x)**: bug fixes, minor GUI changes, minor added functionalities. Re-installation is needed.

Micro changes **(x.x.1)**: external tools update (e.g. Vina, Smina, ADFR). The *Check for Updates* function can be used to update the plugin. See the [User's Guide](https://github.com/paiardin/DockingPie/releases/download/versioning/DockingPie_User_Guide.pdf "DockingPie User's guide") for further details. 
    

## References

[1] DeLano, WL. (2002). “The PyMOL Molecular Graphics System on World Wide Web.” CCP4 Newsletter On Protein Crystallography

[2] Ruiz-Carmona, S., Alvarez-Garcia, D., Foloppe, N., Garmendia-Doval, A. B., Juhos, S., Schmidtke, P., Barril, X., Hubbard, R. E., & Morley, S. D. (2014). rDock: a fast, versatile and open source program for docking ligands to proteins and nucleic acids. PLoS computational biology, 10(4), e1003571. https://doi.org/10.1371/journal.pcbi.1003571

[3] O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461

[4] Koes, David Ryan et al. “Lessons learned in empirical scoring with smina from the CSAR 2011 benchmarking exercise.”Journal of chemical information and modeling vol. 53,8 (2013): 1893-904. doi:10.1021/ci300604z

[5] Ravindranath PA, Forli S, Goodsell DS, Olson AJ, Sanner MF. AutoDockFR: Advances in Protein-Ligand Docking with Explicitly Specified Binding Site Flexibility. PLoS Comput Biol. 2015;11(12):e1004586. Published 2015 Dec 2. doi:10.1371/journal.pcbi.1004586

[6] Morris, G. M., Huey, R., Lindstrom, W., Sanner, M. F., Belew, R. K., Goodsell, D. S., & Olson, A. J. (2009). AutoDock4 and AutoDockTools4: Automated docking with selective receptor flexibility. Journal of Computational Chemistry, 30(16), 2785–2791.

[7] O'Boyle, N.M., Banck, M., James, C.A. et al. Open Babel: An open chemical toolbox. J Cheminform 3, 33 (2011). https://doi.org/10.1186/1758-2946-3-33

[8] Meli, R., Biggin, P.C. spyrmsd: symmetry-corrected RMSD calculations in Python. J Cheminform 12, 49 (2020). https://doi.org/10.1186/s13321-020-00455-2

[9] https://sourceforge.net/projects/sdsorter/

[10] Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. Retrieved from https://docs.anaconda.com/
 
 
## Acknowledgments 

The authors are grateful to Giacomo Janson for assistance with Python.
The authors are grateful to Vedran Miletić for assistance with RxDock.
