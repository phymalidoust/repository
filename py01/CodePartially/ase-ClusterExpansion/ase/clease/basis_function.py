"""Module for setting up pseudospins and basis functions."""
import numpy as np
import math


class BasisFunction(object):
    """Base class for all Basis Functions."""

    def __init__(self, unique_elements):
        self.unique_elements = sorted(unique_elements)
        self.num_unique_elements = len(unique_elements)
        if self.num_unique_elements < 2:
            raise ValueError("Systems must have more than 1 type of element.")
        self.spin_dict = self.get_spin_dict()
        self.basis_functions = self.get_basis_functions()

    def get_spin_dict(self):
        """Get spin dictionary."""
        raise NotImplementedError("get_spin_dict has to be implemented in "
                                  "derived classes!")

    def get_basis_functions(self):
        """Get basis function."""
        raise NotImplementedError("get_basis_functions has to be implemented "
                                  "in derived classes!")

    def customize_full_cluster_name(self, full_cluster_name):
        """Customize the full cluster names. Default is to do nothing."""
        return full_cluster_name


class Sanchez(BasisFunction):
    """Pseudospin and basis function from Sanchez et al.

    Sanchez, J. M., Ducastelle, F., & Gratias, D. (1984).
    Generalized cluster description of multicomponent systems.
    Physica A: Statistical Mechanics and Its Applications, 128(1-2), 334-350.
    """

    def __init__(self, unique_elements):
        BasisFunction.__init__(self, unique_elements)
        if self.num_unique_elements > 6:
            raise ValueError("Only systems consisting up to 6 types of "
                             "elements are currently supported for the scheme "
                             "by Sanchez et al.")

    def get_spin_dict(self):
        """Define pseudospins for all consistuting elements."""
        # Find odd/even
        spin_values = []
        if self.num_unique_elements % 2 == 1:
            highest = (self.num_unique_elements - 1) / 2
        else:
            highest = self.num_unique_elements / 2
        # Assign spin value for each element
        while highest > 0:
            spin_values.append(highest)
            spin_values.append(-highest)
            highest -= 1
        if self.num_unique_elements % 2 == 1:
            spin_values.append(0)

        spin_dict = {}
        for x in range(self.num_unique_elements):
            spin_dict[self.unique_elements[x]] = spin_values[x]
        return spin_dict

    def get_basis_functions(self):
        """Create basis functions to guarantee the orthonormality."""
        # coeff_c = [c_0^1, c_1^1, c_0^2, c_1^2, c_2^2]^T
        coeff_c = np.array([
            [0.0, np.sqrt(2), -5. / 3, -1 * np.sqrt(10. / 7), -np.sqrt(2)],
            [0.0,  -3 / np.sqrt(2), 2. / 3, np.sqrt(5. / 14),
             3. / (7 * np.sqrt(2))],
            [0.0, 0.0, 0.0, 3 * np.sqrt(2. / 7), 9 * np.sqrt(1.5) / 5],
            [0.0, 0.0, 0.0, -155 / (12 * np.sqrt(14)),
             -101 / (28 * np.sqrt(6))],
            [0.0, 0.0, 0.0, 5 * np.sqrt(7. / 2) / 12, 7 / (20 * np.sqrt(6))]])

        # coeff_d = [d_0^0, d_0^1, d_1^1, d_0^2, d_1^2, d_2^2]^T
        coeff_d = np.array([
            [1.0, np.sqrt(3. / 2), np.sqrt(2. / 5), 1. / np.sqrt(2),
             -np.sqrt(3. / 14)],
            [0.0, 0.0,  -17 / (3 * np.sqrt(10)), -17 / (6 * np.sqrt(2)),
             -7. / 6],
            [0.0, 0.0, np.sqrt(2.5) / 3, 5 / (6 * np.sqrt(2)), 1. / 6],
            [0.0, 0.0, 0.0, 0.0, 131. / (15 * np.sqrt(14))],
            [0.0, 0.0, 0.0, 0.0, -7 * np.sqrt(7. / 2) / 12],
            [0.0, 0.0, 0.0, 0.0, np.sqrt(7. / 2) / 20]])

        col = self.num_unique_elements - 2
        bf_list = []

        bf = {}
        for key, value in self.spin_dict.items():
            bf[key] = coeff_d[0][col] * value
        bf_list.append(bf)

        if self.num_unique_elements > 2:
            bf = {}
            for key, value in self.spin_dict.items():
                bf[key] = (coeff_c[0][col] + (coeff_c[1][col] * value**2))
            bf_list.append(bf)

        if self.num_unique_elements > 3:
            bf = {}
            for key, value in self.spin_dict.items():
                bf[key] = coeff_d[1][col] * value
                bf[key] += coeff_d[2][col] * (value**3)
            bf_list.append(bf)

        if self.num_unique_elements > 4:
            bf = {}
            for key, value in self.spin_dict.items():
                bf[key] = coeff_c[2][col] + (coeff_c[3][col] * (value**2))
                bf[key] += coeff_c[4][col] * (value**4)
            bf_list.append(bf)

        if self.num_unique_elements > 5:
            bf = {}
            for key, value in self.spin_dict.items():
                bf[key] = coeff_d[3][col] * value
                bf[key] += coeff_d[4][col] * (value**3)
                bf[key] += coeff_d[5][col] * (value**5)
            bf_list.append(bf)

        return bf_list


class VandeWalle(BasisFunction):
    """Pseudospin and basis function from van de Walle.

    van de Walle, A. (2009).
    Multicomponent multisublattice alloys, nonconfigurational entropy and other
    additions to the Alloy Theoretic Automated Toolkit. Calphad, 33(2),
    266-278.
    """

    def __init__(self, unique_elements):
        BasisFunction.__init__(self, unique_elements)

    def get_spin_dict(self):
        """Define pseudospins for all consistuting elements."""
        spin_values = list(range(self.num_unique_elements))
        spin_dict = {}
        for x in range(self.num_unique_elements):
            spin_dict[self.unique_elements[x]] = spin_values[x]
        return spin_dict

    def get_basis_functions(self):
        """Create basis functions to guarantee the orthonormality."""
        alpha = list(range(1, self.num_unique_elements))
        bf_list = []

        for a in alpha:
            bf = {}
            for key, value in self.spin_dict.items():
                var = 2 * np.pi * math.ceil(a/2.) * value
                var /= self.num_unique_elements
                if a % 2 == 1:
                    bf[key] = -np.cos(var) + 0.
                else:
                    bf[key] = -np.sin(var) + 0.

            # normalize the basis function
            sum = 0
            for key, value in self.spin_dict.items():
                sum += bf[key] * bf[key]
            normalization_factor = np.sqrt(self.num_unique_elements / sum)

            for key, value in bf.items():
                bf[key] = value * normalization_factor

            bf_list.append(bf)

        return bf_list


def kronecker(i, j):
    """Kronecker delta function."""
    if i == j:
        return 1
    return 0


class Sluiter(BasisFunction):
    """Pseudospin and basis function from Sluiter.

    Zhang, Xi, and Marcel HF Sluiter.
    Cluster expansions for thermodynamics and kinetics of multicomponent
    alloys.
    Journal of Phase Equilibria and Diffusion 37.1 (2016): 44-52.
    """

    def __init__(self, unique_elements, reduntant_element="auto"):
        if reduntant_element == "auto":
            self.redundant_element = sorted(unique_elements)[0]
        else:
            self.redundant_element = reduntant_element
        BasisFunction.__init__(self, unique_elements)

    def get_spin_dict(self):
        """Define pseudospins for all consistuting elements."""
        spin_values = list(range(self.num_unique_elements))
        spin_dict = {}
        for x in range(self.num_unique_elements):
            spin_dict[self.unique_elements[x]] = spin_values[x]
        return spin_dict

    def get_basis_functions(self):
        """Create orthonormal basis functions.

        Due to the constraint that any site is occupied by exactly one element,
        we only need to track N-1 species if there are N species.
        Hence, the first element specified is redundant, and will not
        have a basis function.
        """
        bf_list = []
        num_bf = self.num_unique_elements
        for bf_num in range(num_bf):
            if self.unique_elements[bf_num] == self.redundant_element:
                continue
            new_bf = {symb: float(kronecker(i, bf_num))
                      for i, symb in enumerate(self.unique_elements)}
            bf_list.append(new_bf)
        return bf_list

    def _decoration2element(self, dec_num):
        """Get the element with its basis function equal to 1."""
        bf = self.bf_list[dec_num]
        for k, v in bf.items():
            if v == 1:
                return k
        raise ValueError("Did not find any element where the value is 1.")

    def customize_full_cluster_name(self, full_cluster_name):
        """Translate the decoration number to element names."""
        dec = full_cluster_name.rpartition("_")[1]
        name = full_cluster_name.rpartition("_")[0]
        new_dec = ""
        for decnum in dec:
            element = self._decoration2element(int(decnum))
            new_dec += "{}".format(element)
        return name + "_" + new_dec
