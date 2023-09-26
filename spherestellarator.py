# KNOWN ERROR, default sphere cad model has preexisting slice
# Syn Hubbard (CNERG)


# Import statements for calculations and export
import cadquery as cq
import math
from cadquery import exporters
import read_vmec


# Create the spherical surface so it can be split accordingly
def getBody (radius):
    return cq.Workplane().sphere(radius)

# Export
def export(sphericalSurface):
    exporters.export(sphericalSurface, './sphericalSurface.stl')

def main():
    vmec = read_vmec.vmec_data(r'C:\Users\vich2\Downloads\plas_eq.nc')
    # x, y, z =vmec.vmec2xyz(s, theta, phi) -> s = 1
    # stellarator = parametricSurface((x, y, z),20,0.0,1.0,1e-3,1,3,(1,1,1))
    # show_object(stellarator)
    # Define parameters for splitting, in reality, it would be nice for this to be inputted by the user
    slices = 20  # replace integer literals -> prompt for slices
    radius = 1
    shellThickness = 0.05
    thetaIncrement = (2 * math.pi) / slices
    totalHeight = radius * 2  # Total height of the sphere
    heightIncrement = totalHeight / slices
    Tcos = 0.0
    Tsin = 0.0
    sphericalSurface = getBody(radius).faces("+Z").shell(shellThickness) # Using the +Z axis, shell it to make hollow

    # In the longitudinal direction, split sphere
    for index in range(slices):
        if (index == 0):
            sphericalSurface.rotateAboutCenter((2 * math.pi/slices),0,0)
        angle = (index * thetaIncrement)
        Tcos = math.cos(angle)
        Tsin = math.sin(angle)
        splittingPlaneLong = cq.Face.makePlane(None, None, cq.Vector(0, 0, 0), cq.Vector(Tcos, Tsin, 0))
        sphericalSurface = sphericalSurface.split(splittingPlaneLong)

    # In the latitudinal direction, split sphere with evenly spaced slices
    for index in range(slices):  # Start from 1 to avoid splitting at the top and bottom poles
        height = -radius + (index * heightIncrement)  # Adjust height based on the increment
        splittingPlaneLat = cq.Face.makePlane(None, None, cq.Vector(0, 0, height), cq.Vector(0, 0, radius))
        sphericalSurface = sphericalSurface.split(splittingPlaneLat)

    show_object(sphericalSurface)
    export(sphericalSurface)

main()