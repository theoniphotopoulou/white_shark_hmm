# white_shark_hmm
HMM analysis of ARGOS satellite tracking data from 33 white sharks tagged in South Africa

This repository contains the functions and scripts used for analysing white shark tracking data using hidden Markov models. The rationale for the analysis is presented in detail in "Sex and size influence the spatiotemporal distribution of white sharks, with implications for interactions with fisheries and spatial management in the southwest Indian Ocean" by Alison Kock, Amanda T Lombard, Ryan Daly, Victoria Goodall, Michael Meÿer, Ryan Johnson, Chris Fischer, Pieter Koen, Dylan Irion, Enrico Gennari, Alison Towner, Oliver Jewell, Charlene da Silva, Matt Dicken, Malcolm Smale, Theoni Photopoulou (2022) published in the journal Frontiers in Marine Science https://doi.org/10.3389/fmars.2022.811985

The data are read in and cleaned in script "01satellite_data.R", the HMMs are fitted in script "02momentuHMM_models.R" and the plots are produced in script "03ARGOSshark_plot_results.R".

The data, output and figures are not available in this repo while the work is under review for publication. The dataset used in the analysis is available at the following URL <https://zenodo.org/record/5575189> 

Please contact Theoni Photopoulou at <theoni.photopoulou@gmail.com> with any questions about the code and the data format.

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
