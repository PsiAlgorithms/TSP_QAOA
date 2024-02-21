import numpy as np


def path_rev(path, k):
    """
    Reverse a path.

    Parameters
    ----------
    path : str
        Binary sequence of a path.
    k : int
        Number of qubits needed to encode cities.

    Returns
    -------
    path_reversed : str
        Revense of the input path.

    """
    path_reversed = path[-k:]
    for i in range(len(path) // k):
        path_reversed += path[-k * (i + 1) : -k * i]
    return path_reversed


def check_validity(sequence, k, encoding):
    """
    Check if a binary sequence is a valid path.

    Parameters
    ----------
    sequence : str
        A binary sequence.
    k : int
        Number of qubits needed to encode cities.
    encoding : str
        Type of encoding.

    Returns
    -------
    bool
        Weather a sequence is a path or not.

    """
    triseq = []
    for i in range(len(sequence) // k):
        triseq.append(sequence[k * i : k * (i + 1)])
        if encoding == "n2" and triseq[-1].count("1") != 1:
            return False
    return True if len(np.unique(triseq)) == len(sequence) // k and "0" * k not in triseq else False


def euclidean_distance(coord1, coord2):
    return np.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2)


def path_distance(coord, path, k, encoding):
    """
    Calculate the distance of a path given the coordinates.

    Parameters
    ----------
    coord : numpy.ndarray
        set of coordinates.
    path : str
        binary sequence of a path.
    k : int
        Number of qubits needed to encode cities.
    encoding : str
        Type of encoding.

    Returns
    -------
    distance : float
        the distance of the input path.

    """
    path_int = [0]

    if encoding == "n2":
        for i in range(len(path) // k):
            path_int.append(path[k * i : k * (i + 1)][::-1].index("1") + 1)
    else:
        for i in range(len(path) // k):
            path_int.append(int(path[k * i : k * (i + 1)], 2))

    path_int.append(0)
    coord = list(coord)
    coord.append(coord[0])
    distance = 0
    for i in range(len(path_int) - 1):
        distance += euclidean_distance(coord[path_int[i]], coord[path_int[i + 1]])
    return distance


def generate_binary_nums(nqubit):
    binary_numbers = []
    for i in range(2**nqubit):
        binary_str = format(i, f"0{nqubit}b")
        binary_numbers.append(binary_str)
    return binary_numbers
