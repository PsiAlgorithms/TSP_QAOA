# TSP_QAOA

This repository contains the implementation of the Traveling Salesman Problem with Quantum Approximate Optimization Algorithm.

## Repository overview

All neceassry packages and libraries are inclunded in requirements.txt.

main.py contains the code to create TSP_QAOA objects with both n^2 and nlogn encodings, and runs the simulation and result analysis.

We have also included our test coordinates witch can be seen in plot_coordinates.ipynb.

hamiltonian.py contains the problem hamiltonian of both n^2 and nlogn encodings, with coefficient.py contatning the coeffients in both encodings' formulations.

## Cite

The respective paper will soon be published.

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
