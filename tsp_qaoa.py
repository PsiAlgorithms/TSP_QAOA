import numpy          as np
import pickle         as pkl
import scipy.optimize as opt

from qiskit.quantum_info    import Statevector
from qiskit.circuit         import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.compiler        import transpile
from qiskit_aer.backends    import statevector_simulator

from utilities import path_distance, path_rev, check_validity, generate_binary_nums
import coefficient
import hamiltonian


class TSP_QAOA:
    def __init__(self, encoding, coordinates):
        if encoding not in ["n2", "nlogn"]:
            raise NotImplementedError("Sorry, encoding should be 'n2' or 'nlogn'.")

        self.encoding = encoding
        self.coordinates = coordinates
        self.n = coordinates.shape[0]  # Number of cities.

        if encoding == "n2":
            self.k = self.n - 1
            self.nqubit = (self.n - 1) ** 2
            self.coeff_generator = coefficient.coeff_n2
            self.problem_hamlitonian = hamiltonian.hamiltonian_n2
        else:
            self.k = int(np.ceil(np.log2(self.n)))
            self.nqubit = self.k * (self.n - 1)
            self.coeff_generator = coefficient.coeff_nlogn
            self.problem_hamlitonian = hamiltonian.hamiltonian_nlogn

        self.mixing_hamiltonian = hamiltonian.mixer
        self.backend = statevector_simulator.StatevectorSimulator()

    def generate_coeffs(self):
        self.coeff = self.coeff_generator(self.coordinates)

    def generate_hamiltonians(self):
        self.hp = self.problem_hamlitonian(self.coeff, self.nqubit, self.n, self.k)
        self.hm = self.mixing_hamiltonian(self.nqubit)

    def UH(self, gamma):
        # Hamiltonian unitary
        evolution = PauliEvolutionGate(self.hp, gamma)
        evolution.label = r"$\exp(-i\gamma H_{P})$"
        return evolution

    def UM(self, beta):
        # Mixing unitary
        evolution = PauliEvolutionGate(self.hm, beta)
        evolution.label = r"$\exp(-i\beta H_{M})$"
        return evolution

    def create_qaoa_circ(self, theta):
        """
        Creates a parametrized QAOA circuit.

        Parameters
        ----------
        theta : list
            Unitary parameters.

        Returns
        -------
        qc : QuantumCircuit
            Qiskit circuit.

        """
        n_layers = len(theta) // 2  # number of alternating unitaries
        beta = theta[:n_layers]
        gamma = theta[n_layers:]

        qreg = []
        for i in range(self.n - 1):
            qreg.append(QuantumRegister(self.k, f"(T{i+2})"))
        qc = QuantumCircuit(*qreg)

        # initial_state
        qc.h(range(self.nqubit))

        # layers
        for layer_index in range(n_layers):
            # problem unitary
            qc.append(self.UH(gamma[layer_index]), range(self.nqubit))
            # mixing unitary
            qc.append(self.UM(beta[layer_index]), range(self.nqubit))

        return qc

    def execute_circuit(self, theta):
        """
        simulate the circuit.

        Parameters
        ----------
        theta : list
            A list containing beta and gamma to be evaluated by the circuit.

        Returns
        -------
        TYPE
            DESCRIPTION.
        """
        qc = self.create_qaoa_circ(theta)
        qc_transpiled = transpile(qc, self.backend)
        job = self.backend.run(qc_transpiled)
        result = job.result()
        state = result.get_statevector(qc)
        return state

    def objective_layer(self, layer, resx):
        layer = layer
        resx = resx

        def wrapper(theta):
            if layer == 0:
                res_list = theta
            else:
                res_list = np.zeros(2 * (layer + 1))
                res_list[slice(layer, len(res_list), layer + 1)] = theta
                for i in range(layer):
                    res_list[slice(i, len(res_list), layer + 1)] = resx[i]
            return np.real(self.execute_circuit(res_list).expectation_value(self.hp))

        return wrapper

    def fresx(self, layer, resx):
        layer += 1
        res_list = np.zeros(2 * layer)
        for i in range(layer):
            res_list[slice(i, len(res_list), layer)] = resx[i]
        return res_list

    def optimize_objective(
        self,
        nlayer=5,
        theta0=[0, 1],
        bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
        niter=500,
    ):
        """
        Optimize the circuit using basinhopping algorithm.

        Parameters
        ----------
        nlayer : int, optional
            Number of layers. The default is 5.
        theta0 : list, optional
            Starting point for optimization. The default is [0, 1].
        bounds : tuple, optional
            Boundries for beta and gamma. The default is ((0, 2 * np.pi), (0, 2 * np.pi)).
        niter : int, optional
            Number of optimization iterations. The default is 500.

        Returns
        -------
        resx : list
            Pairs of beta and gamma values for each layer.

        """
        self.nlayer = nlayer
        self.binary_numbers = generate_binary_nums(self.nqubit)

        full_states = Statevector.from_label(self.binary_numbers[0])
        for i in range(1, len(self.binary_numbers)):
            full_states += Statevector.from_label(self.binary_numbers[i])
        full_states = (1 / np.sqrt(len(self.binary_numbers))) * full_states

        resx = []
        resfun = [np.real(full_states.expectation_value(self.hp))]
        for layer in range(nlayer):
            print("\t layer " + str(layer + 1) + " started")
            res = opt.basinhopping(
                self.objective_layer(layer, resx),
                theta0,
                minimizer_kwargs={"bounds": bounds},
                niter=niter,
            )
            if res.fun < resfun[layer]:
                resx += [list(res.x)]
                resfun += [res.fun]
            else:
                while res.fun >= resfun[layer]:
                    print("\t false")
                    print(res.x, res.fun)
                    res = opt.basinhopping(
                        self.objective_layer(layer, resx),
                        np.random.rand(2) * (2 * np.pi),
                        minimizer_kwargs={"bounds": bounds},
                        niter=niter,
                    )
                resx += [list(res.x)]
                resfun += [res.fun]
            # print("\t layer " + str(layer + 1) + " finished")
            print("\t results of layer " + str(layer + 1) + ":")
            print("\t", res.x, res.fun, "\n")

        self.resx = resx
        return resx

    def rank(self, path_probs_dict, best_path):
        probs_sort_indeces = np.argsort(list(path_probs_dict.values()))[::-1]
        paths_sorted = np.array(list(path_probs_dict.keys()))[probs_sort_indeces]
        rank_best_path = np.where(paths_sorted == best_path)[0]
        rank_best_path_rev = np.where(paths_sorted == path_rev(best_path, self.k))[0]

        return int(min([rank_best_path, rank_best_path_rev]) + 1)

    def results_from_params(self, filename=None):
        valid_paths = [bn for bn in self.binary_numbers if check_validity(bn, self.k, self.encoding)]
        valid_paths_dist_dict = {}
        for vp in valid_paths:
            valid_paths_dist_dict[vp] = path_distance(self.coordinates, vp, self.k, self.encoding)
        best_path = min(valid_paths_dist_dict, key=valid_paths_dist_dict.get)

        state_res = []
        prob_dict_res = []
        xp = []  # expectation
        ar = []  # approximation ratio
        tp = []  # true percentage
        rn = []  # rank

        for layer in range(self.nlayer):
            state_res.append(self.execute_circuit(self.fresx(layer, self.resx[0 : layer + 1])))
            xp.append(np.real(state_res[layer].expectation_value(self.hp)))
            ar.append(xp[layer] / path_distance(self.coordinates, best_path, self.k, self.encoding))
            prob_dict_res.append(state_res[layer].probabilities_dict())
            tp.append(prob_dict_res[layer][best_path] + prob_dict_res[layer][path_rev(best_path, self.k)])
            rn.append(self.rank(prob_dict_res[layer], best_path))

        self.result = {
            "state_res": state_res,
            "prob_dict_res": prob_dict_res,
            "xp": xp,
            "ar": ar,
            "tp": tp,
            "rn": rn,
        }

        if filename:
            with open(f"results/{filename}", "wb") as file:
                pkl.dump(self.result, file)

    def load_results(self, filename):
        with open(f"results/{filename}", "rb") as file:
            self.result = pkl.load(file)
