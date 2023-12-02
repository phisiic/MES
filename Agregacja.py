from WektorP import MacierzHBC, WektorP
from MacierzH import no_integration_nodes
from tabulate import tabulate


class MacierzHGlobalna:
    def __init__(self, no_elements, no_nodes, h_matrices):
        self.no_nodes = no_nodes
        self.h_matrices = h_matrices
        self.elements = []
        self.element_IDs = [[0] * 4 for _ in range(no_elements)]
        self.h_matrix_global = [[0] * no_nodes for _ in range(no_nodes)]

        for matrix in self.h_matrices:
            self.elements.append(matrix.element)

        for i in range(no_elements):
            for j in range(4):
                self.element_IDs[i][j] = self.elements[i].connected_nodes[j].node_id

        for k in range(no_elements):
            h_matrix = self.h_matrices[k].total_matrix
            for i in range(4):
                for j in range(4):
                    global_row = self.element_IDs[k][i]
                    global_col = self.element_IDs[k][j]
                    self.h_matrix_global[global_row - 1][global_col - 1] += h_matrix[i][j]

    def print_global_matrix(self):
        headers = [""] + list(range(1, self.no_nodes + 1))
        table = [[i + 1] + row for i, row in enumerate(self.h_matrix_global)]
        print(tabulate(table, headers=headers, tablefmt="fancy_grid"))


class WektorPGlobalny:
    def __init__(self, no_elements, no_nodes, p_vectors):
        self.p_vectors = p_vectors
        self.p_vector_global = [0 for _ in range(no_nodes)]

        self.elements = []
        self.element_IDs = [[0] * 4 for _ in range(no_elements)]

        for vector in self.p_vectors:
            self.elements.append(vector.element)

        for i in range(no_elements):
            for j in range(4):
                self.element_IDs[i][j] = self.elements[i].connected_nodes[j].node_id

        for k in range(no_elements):
            p_vector = self.p_vectors[k].p_vector
            for i in range(4):
                global_index = self.element_IDs[k][i] - 1  # Adjusting for 0-based indexing
                self.p_vector_global[global_index] += p_vector[i]

    def print_global_vector(self):
        headers = ["Node ID", "P Vector Value"]
        table = [[i + 1, value] for i, value in enumerate(self.p_vector_global)]
        print(tabulate(table, headers=headers, tablefmt="fancy_grid"))

