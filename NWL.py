import numpy as np
from scipy.optimize import direct
import h5py
import matplotlib.pyplot as plt
import read_vmec


def root_problem(theta, vmec, wall_s, phi, pt):
    """
    """
    # Compute first wall point
    wall_pt = np.array(vmec.vmec2xyz(wall_s, theta, phi))
    # Convert from m to cm
    wall_pt = wall_pt*100
    # Compute L2 norm of difference between original and computed points
    diff = np.linalg.norm(pt - wall_pt)

    return diff


def find_coords(vmec, wall_s, phi, pt):
    """
    """
    # Root-solve for the poloidal angle
    theta = direct(
        root_problem,
        bounds = [(0., 2*np.pi)],
        args = (vmec, wall_s, phi, pt)
    )
    # Extract angle
    theta = theta.x

    return theta


def flux_coords(vmec, wall_s, pt):
    """
    """
    # Extract Cartesian coordinates of original point
    x, y, z = pt
    # Compute toroidal angle
    phi = np.arctan2(y, x)
    # Compute poloidal angle
    theta = find_coords(vmec, wall_s, phi, pt)

    return phi, theta


def extract_coords(source_file):
    """
    """
    # Load source file
    file = h5py.File(source_file, 'r')
    # Extract source information
    dataset = file['source_bank']
    # Extract coordinates of particle crossings
    coords = dataset['r']

    return coords


# Define source file
source_file = '/Users/synh/Desktop/parastell/parametric_stellarator/surface_source.h5'
# Extract coordinates of particle crossings
coords = extract_coords(source_file)
# Define plasma equilibrium VMEC file
plas_eq = '/Users/synh/Desktop/parastell/parametric_stellarator/plas_eq.nc'
# Load plasma equilibrium data
vmec = read_vmec.vmec_data(plas_eq)
# Define closed flux surface extrapolation for first wall
wall_s = 1.2
# Define toroidal extent of one period
tor_ext = 90.
# Define poloidal extent of plasma equilibrium
pol_ext = 360.

# Define number of bins for toroidal and poloidal angles
num_phi = 101
num_theta = 101
# Define toroidal and poloidal angle bins
phi_bins_cent = np.linspace(0.0, tor_ext, num = num_phi)
theta_bins_cent = np.linspace(-pol_ext/2, pol_ext/2, num = num_theta)
"""
*For Syn's reference*

Actual bin boundaries are located at midpoints of bins. E.g.,

phi_bins_cent = [0, 30, 60, 90]

gives bin boundaries

phi_bins_bounds = [0, 15, 45, 75, 90]
"""
# Initialize count matrix
count_mat = np.zeros((num_phi, num_theta))

# Loop over coordinates of surface crossings
for pt in coords:
    # Extract Cartesian coordinates and format as array
    pt = [pt['x'], pt['y'], pt['z']]
    pt = np.array(pt)
    # Compute toroidal and poloidal angles of first wall point
    phi, theta = flux_coords(vmec, wall_s, pt)
    # Convert angles from radians to degrees
    phi = np.rad2deg(phi)
    theta = np.rad2deg(theta)
    # Shift angles to fit in bins
    if theta > pol_ext/2:
        theta = theta - pol_ext

    # Loop over toroidal angle bins
    for i, phi_bin in enumerate(phi_bins_cent):
        # Conditionally...
        if np.abs(phi - phi_bin) <= tor_ext/(num_phi - 1)/2:
            # Loop over poloidal angle bins
            for j, theta_bin in enumerate(theta_bins_cent):
                # Conditionally...
                if np.abs(theta - theta_bin) <= pol_ext/(num_theta - 1)/2:
                    count_mat[i,j] += 1

# Define neutron energy (eV)
n_energy = 14.1e6
# Define eV to joules constant
eV2J = 1.60218e-19
# Compute total neutron source strength (n/s)
SS = 1.873228252667338e+20
# Define joules to megajoules constant
J2MJ = 1e-6
# Define number of source particles
num_parts = len(coords)

# Compute NWL
NWL_mat = count_mat*n_energy*eV2J*SS*J2MJ/num_parts

# Compute surface area of each bin
"""
This is where Syn's work will go to compute surface area of each bin :)
"""
# Adjust bin boundaries so that it accurately describes endpoints
# Initialize arrays with beginning points, then interate over mid points
num_phi_bins = len(phi_bins_cent) - 2
num_theta_bins = len(theta_bins_cent) - 2
surf_area_mat = []
phi_bins_boundaries = [0.0]
theta_bins_boundaries = [-pol_ext/2]
vectors_1 = []
vectors_2= []

# Compute vectors defining the bin in Cartesian coordinates
# TODO: CHECK TO MAKE SURE VECTORS ARE ACCURATE
for i in range(num_phi_bins):
    if i == num_phi_bins - 1:
        phi_bin_boundary = 2*np.pi
    else:
        phi_bin_boundary = (phi_bins_cent[i] + phi_bins_cent[i+1]) / 2
    for j in range(num_theta_bins):
        if j == num_theta_bins - 1:
            theta_bin_boundary = 2*np.pi
        else:
            theta_bin_boundary = (theta_bins_cent[j] + theta_bins_cent[j+1]) / 2
            
    vec1 = np.array([wall_s * np.cos(phi_bin_boundary), wall_s * np.sin(phi_bin_boundary), 0.0])
    vec2 = np.array([-wall_s * np.sin(phi_bin_boundary) * np.sin(theta_bin_boundary),
                         wall_s * np.cos(phi_bin_boundary) * np.sin(theta_bin_boundary),
                         wall_s * np.cos(theta_bin_boundary)])
    vectors_1.append(vec1)
    vectors_2.append(vec2)

for i in range(num_phi_bins):
    area_vector = np.cross(vectors_1[i], vectors_2[i])
    surf_area = np.linalg.norm(area_vector)
    surf_area_mat.append(surf_area)

# Plot NWL normalized by surface area
NWL_mat_normalized = NWL_mat / surf_area_mat

# Plot NWL
levels = np.linspace(np.min(NWL_mat), np.max(NWL_mat), num = 101)
fig, ax = plt.subplots()
CF = ax.contourf(phi_bins_cent, theta_bins_cent, NWL_mat, levels = levels)
cbar = plt.colorbar(CF)
cbar.ax.set_ylabel('NWL (MW)')
plt.xlabel('Toroidal Angle (degrees)')
plt.ylabel('Poloidal Angle (degrees)')
fig.savefig('NWL.png')
