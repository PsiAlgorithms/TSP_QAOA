# TSP_QAOA

This repository implements the Traveling Salesman Problem (TSP) using the Quantum Approximate Optimization Algorithm (QAOA). QAOA is a variational algorithm for solving optimization problems on quantum computers, and this project explores its application to finding optimal routes for the TSP.

## The Traveling Salesman Problem (TSP)

The Traveling Salesman Problem is a classic problem in combinatorial optimization. It asks the following question: given a list of cities and the distances between them, what is the shortest possible route that visits each city exactly once and returns to the starting city? The TSP has applications in various fields including logistics, manufacturing, and genetic sequencing. Due to its simple formulation but complex nature, finding the optimal solution to the TSP is an NP-hard problem, meaning that the time required to solve it grows exponentially with the number of cities.

## Repository Structure

**Requirements:**

* All necessary packages and libraries are listed in `requirements.txt`. Install them using `pip install -r requirements.txt`.

**Files:**

<p align="center">
  <img alt="Repository Structure" src="./images/repo_structure.svg">
</p>

* `main.ipynb`: Includes the example usage. it creates TSP_QAOA objects using both $n^2$ and $n \log n$ encodings, runs the QAOA simulation, and analyzes the results.
* `tsp_qaoa.py`: Contains the core implementation of the project through the `TSP_QAOA` class. This class allows users to solve the TSP using QAOA with two different encodings.
* `hamiltonian.py`: Defines the problem Hamiltonian for both encodings.
* `coefficient.py`: Calculates the coefficients used in the Hamiltonian formulations.
* `plot_coordinates.ipynb`: This Jupyter notebook visualizes our test coordinates in `data/coordinates.pkl`.

## Citation

Cite as:
```
@article{ramezani2024reducing,
  title={Reducing the Number of Qubits from $ n\^{} 2$ to $ n$\backslash$log\_ $\{$2$\}$(n) $ to Solve the Traveling Salesman Problem with Quantum Computers: A Proposal for Demonstrating Quantum Supremacy in the NISQ Era},
  author={Ramezani, Mehdi and Salami, Sadegh and Shokhmkar, Mehdi and Moradi, Morteza and Bahrampour, Alireza},
  journal={arXiv preprint arXiv:2402.18530},
  year={2024}
}
```

## Contributing

This project is under active development. We welcome contributions to this project! Please follow these guidelines:

* Fork the repository and create a pull request for your changes.
* Ensure your code is well-documented and adheres to existing coding style guidelines.
* Include descriptions of your modifications and any relevant test cases.

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
