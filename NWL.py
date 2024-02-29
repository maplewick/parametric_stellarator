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
# Adjust bin boundaries to accurately describe endpoints
# Bin boundaries should have bounds of [0, ....., 90], [-180, ......, 180]
phi_bins_boundaries = np.linspace(0.0, tor_ext, num_phi + 1)
theta_bins_boundaries = np.linspace(-pol_ext/2, pol_ext/2, num_theta + 1)

# Initialize surface area matrix so that is 
surf_area_mat = np.zeros((num_phi, num_theta))

# Compute surface area for each bin
for i in range(num_phi):
    for j in range(num_theta):
        # Convert VMEC coordinates to Cartesian coordinates for the four corners of the bin (i.e. r1, r2, r3, r4)
        # Lower left corner
        pt1 = np.array(vmec.vmec2xyz(wall_s, theta_bins_boundaries[j], phi_bins_boundaries[i]))
        # Lower right corner
        pt2 = np.array(vmec.vmec2xyz(wall_s, theta_bins_boundaries[j], phi_bins_boundaries[i + 1]))
        # Upper left corner
        pt3 = np.array(vmec.vmec2xyz(wall_s, theta_bins_boundaries[j + 1], phi_bins_boundaries[i]))
        # Upper right corner
        pt4 = np.array(vmec.vmec2xyz(wall_s, theta_bins_boundaries[j + 1], phi_bins_boundaries[i + 1]))

        # Compute vectors along the bin edges
        vector1 = pt4 - pt3
        vector2 = pt1 - pt3
        vector3 = pt4 - pt2
        vector4 = pt1 - pt2

        # Second set of vectors used for aggregation 
        vector5 = pt3 - pt4
        vector6 = pt2 - pt4
        vector7 = pt2 - pt1
        vector8 = pt3 - pt1

        # Compute the cross product of the vectors to find the area of cross section (set 1)
        cross_prod1 = np.cross(vector1, vector2)
        cross_prod2 = np.cross(vector3, vector4)
        cross_prod3 = np.cross(vector5, vector6)
        cross_prod4 = np.cross(vector7, vector8)

        # Since each bin has a triangular area it is one half the bin area 
        bin_area1 = 0.5 * np.linalg.norm(cross_prod1)
        bin_area2 = 0.5 * np.linalg.norm(cross_prod2)
        bin_area3 = 0.5 * np.linalg.norm(cross_prod3)
        bin_area4 = 0.5 * np.linalg.norm(cross_prod4)
        bin_area_avg = (bin_area1 + bin_area2 + bin_area3 + bin_area4) / 2

        # Store the computed area in the surface area matrix
        surf_area_mat[i, j] = NWL_mat[i,j] / (bin_area_avg)
        
# Plot Normalized NWL
levels = np.linspace(np.min(surf_area_mat), np.max(surf_area_mat), num = 101)
fig, ax = plt.subplots()
CF = ax.contourf(phi_bins_cent, theta_bins_cent, surf_area_mat, levels = levels)
cbar = plt.colorbar(CF)
cbar.ax.set_ylabel('NWL (MW/m^2)')
plt.xlabel('Toroidal Angle (degrees)')
plt.ylabel('Poloidal Angle (degrees)')
fig.savefig('NWL.png')
