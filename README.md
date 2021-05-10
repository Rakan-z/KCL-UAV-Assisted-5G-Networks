# BEng-Final Year Individual Project
This repository consists of files which include the thesis report and MATLAB programs for my BEng Final Year Individual Project. Each file is described in the below sections.
- Graded 71/100 (Awarded Distinction). 
- Topic description: Unmanned Aerial Vehicles (UAVs), also known as drones, are any aircraft which do not have a human pilot aboard. The use of Unmanned Aerial Vehicle (UAV) assisted 5G communications in the form of flying base stations (FBSs) has the potential to augment terrestrial 5G network operations in a plethora of different use cases such as information dissemination, data collection for Internet of Things (IoTs) and machine-type communications. In such type of use case scenarios, UAV assisted that 5G networks can increase network efficiency and flexibility by providing seamless connectivity at any location (well beyond the coverage area of a single Base Station (BS)) and in that respect, it can ease the growth in high data demand. We can hence envision a scenario where FBS mounted with wireless communication equipment in stationary 5G BS and it can complement and/or enhance the existing cellular network infrastructure. 
- Focus of Individual Project: To optimize the routes of the FBSs that operate on the boundary of two MBSs (Macro Base Stations).
- Proposed two clustering algorithms which can be used to cluster GNs (ground nodes\users) at the edge of a cell (boundary between two MBSs) based on the 2-opt algorithm, planning the routes of the FBSs. Simulated a multi-cell environment using the Monte Carlo method with modelled flashcrowds to investigate results.
- Project Demonstration (powerpoint presentation): https://emckclac-my.sharepoint.com/:p:/g/personal/k1764173_kcl_ac_uk/ETzzTfrfVnRHi0qzlL8oQd8BecgzOIAsE5LbheF1HrNGiQ?e=rb4JbB

Thesis.pdf
-
The project report (thesis paper).

FBS_Boundary_search_lin.m
-
The program which includes the algorithm that linearly clusters the GNs.

FBS_Boundary_search_pw.m
-
The program which includes the algorithm that clusters the GNs using a Piecewise defined function.

READMEthesisFinal.md
-
A guide to the project files.

Simulation_linear_crowd.m
-
The program which simulates a multi-cell environment with a flashcrowd and uses the FBS_Boundary_search_lin algorithm to plan the FBS routes.

Simulation_linear.m
-
The program which simulates a multi-cell environment and uses the FBS_Boundary_search_lin algorithm to plan the FBS routes.

Simulation_piecewise_crowd.m
-
The program which simulates a multi-cell environment with a flashcrowd and uses the FBS_Boundary_search_pw algorithm to plan the FBS routes.

Simulation_piecewise.m
-
The program which simulates a multi-cell environment and uses the FBS_Boundary_search_pw algorithm to plan the FBS routes.

tspsearch.m
-
The 2-opt algorithm used to build the clustering algorithms.

