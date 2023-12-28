from Classes import Node, Element, Global, Grid
from MacierzH import MatrixH, no_integration_nodes
from WektorP import MacierzHBC, WektorP
from Agregacja import MacierzHGlobalna, WektorPGlobalny, MacierzCGlobalna, sum_matrices, sum_vectors
from GaussianElimination import gaussian_elimination
from MacierzC import MacierzC

data = {}
nodes = {}
elements = []
node_section = False
element_section = False
bc_section = False
h_matrices = []
bc_nodes = set()
hbc_matrices = []
summed_matrices = []
p_vectors = []
c_matrices = []

plik = "Test1_4_4.txt"
#plik = "Test2_4_4_MixGrid.txt"

with open(plik, "r") as file:
    for line in file:
        line = line.strip()

        if line.startswith("*Node"):
            node_section = True
            element_section = False
            bc_section = False
            continue

        elif line.startswith("*Element"):
            element_section = True
            node_section = False
            bc_section = False
            continue

        elif line.startswith("*"):
            node_section = False
            element_section = False

        if node_section:
            parts = line.split(',')
            if len(parts) == 3:
                node_id, x, y = [int(parts[0]), float(parts[1]), float(parts[2])]
                node = Node(node_id, x, y)
                nodes[node_id] = node

        if element_section:
            if line:
                element_id, *element_nodes = map(int, line.split(','))
                element = Element(element_id)

                # Find Node objects corresponding to the IDs
                nodes_in_element = [nodes[node_id] for node_id in element_nodes]

                # Add Node objects to the element
                for node in nodes_in_element:
                    element.addNode(node)

                elements.append(element)

        parts = line.split()
        if len(parts) < 3:
            key = parts[0]
            value = parts[1] if len(parts) > 1 else ""
            data[key] = value
        else:
            key = parts[0] + " " + parts[1]
            value = parts[2]
            data[key] = value

# Process BC section after creating Node objects
with open(plik, "r") as file:
    for line in file:
        line = line.strip()

        if line.startswith("*BC"):
            bc_section = True
            node_section = False
            element_section = False
            continue

        if bc_section:
            if line:
                bc_nodes.update(map(int, line.split(',')))

# Modify existing Node objects based on BC information
for node_id in bc_nodes:
    if node_id in nodes:
        nodes[node_id].BC = 1

global_data = Global(
    simTime=int(data.get('SimulationTime', 0)),
    simStepTime=int(data.get('SimulationStepTime', 0)),
    conductivity=int(data.get('Conductivity', 0)),
    alfa=int(data.get('Alfa', 0)),
    tot=int(data.get('Tot', 0)),
    initialTemp=int(data.get('InitialTemp', 0)),
    density=int(data.get('Density', 0)),
    specificHeat=int(data.get('SpecificHeat', 0)),
    nodesNo=int(data.get('Nodes number', 0)),
    elementsNo=int(data.get('Elements number', 0))
)
global_data.print_values()

grid = Grid(global_data.nodesNo, global_data.elementsNo)
for node in nodes.values():
    grid.addNode(node)

for element in elements:
    #element.printElement()
    grid.addElement(element)

    temp_h = MatrixH(element, no_integration_nodes, global_data.conductivity)
    #temp_h.print_total_matrix()
    h_matrices.append(temp_h)

    temp_hbc = MacierzHBC(element, no_integration_nodes, global_data.alfa)
    hbc_matrices.append(temp_hbc.hbc_matrix)

    temp_p_vector = WektorP(element, no_integration_nodes, global_data.alfa, global_data.tot)
    p_vectors.append(temp_p_vector)

    temp_c = MacierzC(global_data.specificHeat, global_data.density, element)
    c_matrices.append(temp_c)

for h_matrix, hbc_matrix in zip(h_matrices, hbc_matrices):
    # Create a new MatrixH instance
    summed_matrix = MatrixH(h_matrix.element, no_integration_nodes, global_data.conductivity)

    # Use the add_hbc_matrix method to add the hbc_matrix to the total_matrix
    summed_matrix.add_hbc_matrix(hbc_matrix)

    # Append the summed matrix to the list
    summed_matrices.append(summed_matrix)

#grid.printGrid()

# count: int = 1
# for matrix in h_matrices:
#     print(f"Matrix H {count}")
#     matrix.print_total_matrix()
#     print()
#     count += 1

# count: int = 1
# for matrix in hbc_matrices:
#     print(f"Matrix HBC {count}")
#     for i in range(4):
#         for j in range(4):
#             print("{:.10f}".format(matrix[i][j]), end="\t")
#         print()
#     count += 1

# count: int = 1
# for vector in p_vectors:
#     print(f"P Vector {count}")
#     vector.print_total_vector()
#     print()
#     count += 1

# # Iterate through the matrices in h_matrices and hbc_matrices
# for h_matrix, hbc_matrix in zip(h_matrices, hbc_matrices):
#     # Initialize a new matrix for the sum
#     summed_matrix = [[0.0] * 4 for _ in range(4)]
#
#     # Add corresponding elements of h_matrix and hbc_matrix
#     for i in range(4):
#         for j in range(4):
#             summed_matrix[i][j] = h_matrix.total_matrix[i][j] + hbc_matrix[i][j]
#
#     # Append the summed matrix to the list
#     summed_matrices.append(summed_matrix)

#print("\n Macierz H Globalna")
h_glob = MacierzHGlobalna(global_data.elementsNo, global_data.nodesNo, summed_matrices)
#macierz_h_globalna.print_global_matrix()

#print("\n Wektor P Globalny")
p_glob = WektorPGlobalny(global_data.elementsNo, global_data.nodesNo, p_vectors)
#wektor_p_globalny.print_global_vector()

#print("\n Macierz C Globalna")
c_glob = MacierzCGlobalna(global_data.elementsNo, global_data.nodesNo, c_matrices)
#macierz_c_globalna.print_global_matrix()
c_glob.divide_matrix_by_dtau(global_data.simStepTime)
#macierz_c_globalna.print_global_matrix()

t0_vector = [global_data.initialTemp] * global_data.nodesNo
current_time = 50
while current_time <= global_data.simTime:
    macierz_h_globalna = h_glob
    wektor_p_globalny = p_glob
    macierz_c_globalna = c_glob

    matrix_c_h_summed = sum_matrices(macierz_c_globalna.c_matrix_global, macierz_h_globalna.h_matrix_global)
    matrix_c_multiplied = macierz_c_globalna.multiply_matrix_by_vector(t0_vector)
    vectors_summed = sum_vectors(matrix_c_multiplied, wektor_p_globalny.p_vector_global)

    solution = gaussian_elimination(matrix_c_h_summed, vectors_summed)
    print(f"Solution in time {current_time}:")
    for value in solution:
        print(f"{value:.2f}", end=", ")
    print()

    t0_vector = solution
    current_time += global_data.simStepTime