import itertools
from qiskit.quantum_info import SparsePauliOp
from scipy.special import binom


def hamiltonian_n2(coeff, nqubit, *args):
    identity_list = ["I"] * nqubit
    ## Hamiltonian of Time Steps 1 and 2
    pauli_list = identity_list.copy()
    pauli_str = "".join(map(str, pauli_list))
    h_opt = SparsePauliOp(pauli_str, coeff[0][1])
    for ncoupling in range(1, 3):
        for p in range(int(binom(nqubit, ncoupling))):
            pauli_list = identity_list.copy()
            for q in range(ncoupling):
                qubit_number = coeff[ncoupling][p][q]
                pauli_list[qubit_number] = "Z"
            pauli_str = "".join(map(str, pauli_list))
            pauli_str = pauli_str[::-1]
            h_opt += SparsePauliOp(pauli_str, coeff[ncoupling][p][-1])
    return h_opt


def hamiltonian_nlogn(coeff, nqubit, n, k):
    def index_list(n):
        return list(map(list, itertools.combinations(range(n), 2)))

    coeff0, coeff1, coeff2 = coeff[0], coeff[1], coeff[2]

    identity_list = ["I"] * nqubit
    ## Hamiltonian of Time Steps 1 and 2
    pauli_list = identity_list.copy()
    pauli_str = "".join(map(str, pauli_list))
    h_opt = SparsePauliOp(pauli_str, coeff0[0][1])
    for ncoupling in range(1, k + 1):
        for p in range(int(binom(k, ncoupling))):
            pauli_list = identity_list.copy()
            for q in range(ncoupling):
                qubit_number = coeff0[ncoupling][p][q]
                pauli_list[qubit_number] = "Z"
            pauli_str = "".join(map(str, pauli_list))
            pauli_str = pauli_str[::-1]
            h_opt += SparsePauliOp(pauli_str, coeff0[ncoupling][p][-1])
    ## Hamiltonian of Time Steps i and j
    for i, j in index_list(n - 1):
        coeff = coeff1 if j == i + 1 else coeff2
        pauli_list = identity_list.copy()
        pauli_str = "".join(map(str, pauli_list))
        h_opt += SparsePauliOp(pauli_str, coeff[0][1])
        for ncoupling in range(1, 2 * k + 1):
            for p in range(int(binom(2 * k, ncoupling))):
                pauli_list = identity_list.copy()
                for q in range(ncoupling):
                    quo = coeff[ncoupling][p][q] // k
                    rem = coeff[ncoupling][p][q] % k
                    if quo == 0:
                        qubit_number = (i * k) + rem
                    else:
                        qubit_number = (j * k) + rem
                    pauli_list[qubit_number] = "Z"
                pauli_str = "".join(map(str, pauli_list))
                pauli_str = pauli_str[::-1]
                h_opt += SparsePauliOp(pauli_str, coeff[ncoupling][p][-1])
    ## Hamiltonian of Time Steps n and 1
    pauli_list = identity_list.copy()
    pauli_str = "".join(map(str, pauli_list))
    h_opt += SparsePauliOp(pauli_str, coeff0[0][1])
    for ncoupling in range(1, k + 1):
        for p in range(int(binom(k, ncoupling))):
            pauli_list = identity_list.copy()
            for q in range(ncoupling):
                qubit_number = k * (n - 2) + coeff0[ncoupling][p][q]
                pauli_list[qubit_number] = "Z"
            pauli_str = "".join(map(str, pauli_list))
            pauli_str = pauli_str[::-1]
            h_opt += SparsePauliOp(pauli_str, coeff0[ncoupling][p][-1])
    return h_opt


def mixer(nqubit):
    identity_list = ["I"] * nqubit
    h_opt = SparsePauliOp("".join(map(str, identity_list)), 0)
    for qubit in range(nqubit):
        pauli_list = identity_list.copy()
        pauli_list[qubit] = "X"
        pauli_str = "".join(map(str, pauli_list))
        pauli_str = pauli_str[::-1]
        h_opt += SparsePauliOp(pauli_str)
    return h_opt
