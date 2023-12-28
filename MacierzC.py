from WektorP import N1, N2, N3, N4
from MacierzH import no_integration_nodes, JacobianMatrix
from ElementUniwersalny import UniversalElement
from Classes import Node, Element
from tabulate import tabulate

el = UniversalElement(no_integration_nodes)

n1 = Node(1, 0, 0, 1)
n2 = Node(2, 0.025, 0, 1)
n3 = Node(3, 0.025, 0.025, 1)
n4 = Node(4, 0, 0.025, 1)

# n1 = Node(1, 0.1, 0.005, 1)
# n2 = Node(2, 0.0546918, 0.005, 1)
# n3 = Node(6, 0.0623899, -0.0326101, 0)
# n4 = Node(5, 0.1, -0.0403082, 1)
elem = Element(1)
elem.addNode(n1)
elem.addNode(n2)
elem.addNode(n3)
elem.addNode(n4)


class MacierzC:
    def __init__(self, specific_heat, density, element):

        self.n_functions = [[0] * 4 for _ in range(no_integration_nodes ** 2)]
        self.c_matrices = []
        self.element = element

        for i in range(no_integration_nodes ** 2):
            self.n_functions[i][0] = N1(el.integration_points[i].x, el.integration_points[i].y)
            self.n_functions[i][1] = N2(el.integration_points[i].x, el.integration_points[i].y)
            self.n_functions[i][2] = N3(el.integration_points[i].x, el.integration_points[i].y)
            self.n_functions[i][3] = N4(el.integration_points[i].x, el.integration_points[i].y)

        self.jacobian_determinants = []
        # Loop over each integration point to calculate and store the determinant
        for i in range(no_integration_nodes ** 2):
            temporary_jacobian = JacobianMatrix(element, no_integration_nodes, i)
            self.jacobian_determinants.append(temporary_jacobian.detJ)

        self.weights = []
        for weight in range(len(el.weights)):
            tmp = el.weights[weight]
            self.weights.append(tmp)

        for i in range(no_integration_nodes ** 2):
            matrix_c_temp = [[0] * 4 for _ in range(4)]

            for row in range(4):
                for col in range(4):
                    matrix_c_temp[row][col] += (self.n_functions[i][row] * self.n_functions[i][col] * specific_heat * density * self.jacobian_determinants[i])

            self.c_matrices.append(matrix_c_temp)

        self.total_matrix = self.sum_matrices()

    def sum_matrices(self):
        total_matrix = [[0] * 4 for _ in range(4)]

        for i, matrix in enumerate(self.c_matrices):
            for row in range(4):
                for col in range(4):
                    total_matrix[row][col] += (matrix[row][col] * self.weights[i].x * self.weights[i].y)

        return total_matrix


    def print_N_functions(self):
        for i in range(no_integration_nodes ** 2):
            for j in range(4):
                value = self.n_functions[i][j]
                print(f"{value: .6f}", end="\t")
            print()

    def print_c_matrices(self):
        for i, matrix in enumerate(self.c_matrices):
            print(f"Integration Point {i + 1} Matrix:")
            print(tabulate(matrix, tablefmt="fancy_grid", floatfmt=".6f"))
            print()

    def print_total_matrix(self):
        print("Total Matrix C:")
        print(tabulate(self.total_matrix, tablefmt="fancy_grid", floatfmt=".6f"))


# a = MacierzC(700, 7800, elem)
# a.print_N_functions()
# a.print_c_matrices()
# a.print_total_matrix()
