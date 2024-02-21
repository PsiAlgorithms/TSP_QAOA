import numpy as np
import sympy as sp
import itertools
from utilities import euclidean_distance


def distance_matrix(coordinates):
    W = np.array([[euclidean_distance(coord1, coord2) for coord1 in coordinates] for coord2 in coordinates])
    return W


def coeff_n2(coordinates):
    n = coordinates.shape[0]

    def ind2vec(i, j, cr=n):
        if i >= cr or j >= cr:
            raise Exception("index greater than cols/rows")
        return j * cr + i

    W = distance_matrix(coordinates)

    X = np.empty((n, n), dtype=object)
    for i in range(1, n):
        for j in range(1, n):
            X[i, j] = sp.Symbol(f"x_{i}{j}", bool=True)
    X[0, :] = 0
    X[:, 0] = 0
    X[0, 0] = 1

    Y = np.empty(((n - 1) ** 2), dtype=object)
    for i in range(len(Y)):
        Y[i] = sp.Symbol(f"y_{i}")

    Z = np.empty_like(Y)
    for i in range(len(Z)):
        Z[i] = sp.Symbol(f"z_{i}")

    # Hamiltonian
    A = 10 * np.max(W)
    B = 1

    term1 = sp.Integer(0)
    for i in range(n):
        term1 += sp.expand((1 - np.sum(X[i, :])) ** 2)
    term1 = A * term1

    term2 = sp.Integer(0)
    for j in range(n):
        term2 += sp.expand((1 - np.sum(X[:, j])) ** 2)
    term2 = A * term2

    term3 = sp.Integer(0)
    for u in range(n):
        for v in range(n):
            term3_partial = 0
            for i in range(n - 1):
                term3_partial += X[u, i] * X[v, i + 1]
            term3_partial += X[u, n - 1] * X[v, 0]
            term3 += W[u, v] * term3_partial
    term3 = B * term3

    H = sp.expand(term1 + term2 + term3)

    for xr in X[1:, 1:]:
        for x in xr:
            H = H.subs(x**2, x)

    # Function of Y
    for i, xr in enumerate(X[1:, 1:]):
        for j, x in enumerate(xr):
            H = H.subs(x, Y[ind2vec(i, j, cr=n - 1)])

    # Function of Z
    for i, y in enumerate(Y):
        H = H.subs(y, 1 / 2 * (1 - Z[i]))
    H = H.expand()

    # coeffs
    coeff = []

    coeff.append([0, sp.Poly(H).coeffs()[-1]])

    coeff.append([])
    for i in range(len(Z)):
        coeff_z = sp.Poly(H.coeff(Z[i])).coeffs()[-1]
        coeff[-1].append([i, coeff_z])

    coeff.append([])
    for i in range(len(Z)):
        for j in range(i + 1, len(Z)):
            coeff[-1].append([i, j, H.coeff(Z[i] * Z[j])])

    return coeff


def coeff_nlogn(coordinates):
    n = coordinates.shape[0]
    k = int(np.ceil(np.log2(n)))

    def generate_coupling(m):
        coupling = ["0"]
        coupling.extend(map(str, range(0, m)))
        for k in range(2, n + 1):
            comb = itertools.combinations(coupling[1 : m + 1], k)
            for c in comb:
                coupling.append("".join(c))
        return coupling

    def generate_inputs(m):
        binary_numbers = list(itertools.product([0, 1], repeat=m))
        res = [x + y for x, y in itertools.product(binary_numbers, repeat=2)]
        for i in range(len(res)):
            res[i] = res[i][::-1]
        return res

    def coupling_coefficient(inpt, cpl):
        result = [1]
        for s in cpl[1:]:
            x = 1
            for digit in s:
                x *= inpt[int(digit)]
            result.append(x)
        return result

    def ising(a):
        return [1 if val == 0 else -1 for val in a]

    def generate_coeff(m, cpl, a):
        coeff = [[0, a[0]]]
        for k in range(1, m + 1):
            t1 = []
            for i, s in enumerate(cpl[1:]):
                if len(s) != k:
                    continue
                t2 = []
                for digit in s:
                    t2.append(int(digit))
                t2.append(a[i + 1])
                t1.append(t2)
            coeff.append(t1)
        return coeff

    def coeff_from_coord(coordinates, cpl0, cpl1, Hinv):
        # create distance from coordinates
        d = distance_matrix(coordinates).reshape(-1)

        # penalties
        penalty = 10 * max(d)
        # set penalties
        d1 = d.copy()
        d1[d == 0] = penalty
        d2 = d1.copy()
        d2[d != 0] = 0
        # find coefficients
        a1 = Hinv @ d1
        a2 = Hinv @ d2
        # Hamiltonian neighbors
        Hngb = []
        Hngb.append(a1[0])
        for i, s in enumerate(cpl1[1:]):
            Hngb.append(a1[i + 1])
            for digit in s:
                Hngb[-1] *= sp.Symbol(f"x_{digit}")
        Hngb = np.sum(Hngb)
        # Hamiltonian End Point
        Hep = Hngb
        for i in range(k):
            Hep = Hep.subs(f"x_{i+k}", 1)
        Hep = sp.Poly(Hep)
        # a0
        a0 = [Hep.coeffs()[-1]]
        for i in generate_coupling(k)[1:]:
            x = 1
            for j in i:
                x *= sp.Symbol(f"x_{j}")
            a0.append(Hep.coeff_monomial(x))
        # coefficient
        coeff = [
            generate_coeff(1 * k, cpl0, a0),
            generate_coeff(2 * k, cpl1, a1),
            generate_coeff(2 * k, cpl1, a2),
        ]

        return coeff

    cpl0 = generate_coupling(k)
    cpl1 = generate_coupling(2 * k)

    H = np.zeros((n**2, n**2))
    for i, inpt in enumerate(generate_inputs(k)):
        H[i, :] = coupling_coefficient(ising(inpt), cpl1)
    Hinv = np.linalg.inv(H)

    coeff = coeff_from_coord(coordinates, cpl0, cpl1, Hinv)
    return coeff
