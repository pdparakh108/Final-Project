import numpy as np
import timeit
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define a dot product function used for the rotation operation
def v_dot(a):
    return lambda b: np.dot(a, b)

class LatticeSAW4D:
    def __init__(self, N, l0):
        self.N = N
        self.l0 = l0
        # Initial configuration: a straight line in 4D
        self.init_state = np.dstack((np.arange(N), np.zeros(N), np.zeros(N), np.zeros(N)))[0]
        self.state = self.init_state.copy()
        
        # Define 4D rotation matrices for 90°, 180°, and 270° rotations
        self.rotation_matrices = [
            # 90° rotations
            np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),  # xy-plane 90°
            np.array([[0, 0, -1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),  # xz-plane 90°
            np.array([[0, 0, 0, -1], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]]),  # xw-plane 90°
            np.array([[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),  # yz-plane 90°
            np.array([[1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0]]),  # yw-plane 90°
            np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]]),  # zw-plane 90°

            # 180° rotations
            np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),  # xy-plane 180°
            np.array([[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]),  # xz-plane 180°
            np.array([[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]),  # xw-plane 180°
            np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]),  # yz-plane 180°
            np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]),  # yw-plane 180°
            np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]),  # zw-plane 180°

            # 270° rotations
            np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),  # xy-plane 270°
            np.array([[0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1]]),  # xz-plane 270°
            np.array([[0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [-1, 0, 0, 0]]),  # xw-plane 270°
            np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]]),  # yz-plane 270°
            np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0]]),  # yw-plane 270°
            np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])   # zw-plane 270°
        ]
        
        self.end_to_end_distances = []
        self.radii_of_gyration = []
        self.mean_square_monomer_distances = []
        
    # Define pivot algorithm process where t is the number of successful steps
    def walk(self, t):
        acpt = 0
        # While loop until the number of successful steps reaches t
        while acpt < t:
            pick_pivot = np.random.randint(1, self.N - 1)  # Pick a pivot site
            pick_side = np.random.choice([-1, 1])  # Pick a side
            
            if pick_side == 1:
                old_chain = self.state[0:pick_pivot + 1]
                temp_chain = self.state[pick_pivot + 1:]
            else:
                old_chain = self.state[pick_pivot:]
                temp_chain = self.state[0:pick_pivot]
            
            # Pick a rotation matrix
            rotation_matrix = self.rotation_matrices[np.random.randint(len(self.rotation_matrices))]
            # New chain after rotation
            new_chain = np.apply_along_axis(v_dot(rotation_matrix), 1, temp_chain - self.state[pick_pivot]) + self.state[pick_pivot]

            # Use cdist function of scipy package to calculate the pair-pair distance between old_chain and new_chain
            overlap = cdist(new_chain, old_chain)
            overlap = overlap.flatten()
            
            # Determine whether the new state is accepted or rejected
            if len(np.nonzero(overlap)[0]) != len(overlap):
                continue
            else:
                if pick_side == 1:
                    self.state = np.concatenate((old_chain, new_chain), axis=0)
                elif pick_side == -1:
                    self.state = np.concatenate((new_chain, old_chain), axis=0)
                acpt += 1
        
        # Place the center of mass of the chain on the origin
        self.state = self.l0 * (self.state - np.mean(self.state, axis=0))
        
        # Store end-to-end distance and radius of gyration for calculation of critical exponents
        end_to_end_distance = np.linalg.norm(self.state[-1] - self.state[0])
        self.end_to_end_distances.append(end_to_end_distance**2)
        
        radius_of_gyration = np.mean([np.linalg.norm(position - np.mean(self.state, axis=0))**2 for position in self.state])
        self.radii_of_gyration.append(radius_of_gyration)
        
        # Calculate mean-square distance of a monomer from the end-points
        monomer_distances = np.mean([0.5 * ((np.linalg.norm(self.state[0] - position)**2 + np.linalg.norm(self.state[-1] - position)**2)) for position in self.state])
        self.mean_square_monomer_distances.append(monomer_distances)

# Parameters
l0 = 1    # Bond length (step length)
t = 1000  # Number of pivot steps to perform
N_values = [5, 10, 20, 35, 50, 60, 80, 100, 150, 200, 300, 400, 500, 1000, 2000, 5000, 10000]
num_iter = 1  # Number of iterations for each N value
output_file = "simulation_results_4d.csv"

# Write the headers to the file
with open(output_file, 'w') as f:
    f.write("N,iteration,R_squared,gyration_squared,monomer_distance_squared\n")

for N in N_values:
    for i in range(num_iter):
        chain = LatticeSAW4D(N, l0)
        chain.walk(t)
        
        # Compute mean squared values
        R_squared = chain.end_to_end_distances[0]
        gyration_squared = chain.radii_of_gyration[0]
        monomer_distance_squared = chain.mean_square_monomer_distances[0]
        
        # Save the value to the file
        with open(output_file, 'a') as f:
            f.write(f"{N},{i},{R_squared},{gyration_squared},{monomer_distance_squared}\n")