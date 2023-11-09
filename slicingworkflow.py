# Import statements for calculations and export
import cadquery as cq
import math
from cadquery import exporters as ex
import read_vmec as rv
import OCP
import numpy as np
vmec = rv.vmec_data(r'C:\Users\vich2\Downloads\plas_eq.nc') # Global scope for data

# Export function (exports model as a file)
def export(surface):
    ex.export(surface, './stellarator.step')

# Using the data, the data can be mapped to u and v
# So that the stellarator model can be made parametrically
# with .parametricSurface
# Returns - cartesian coordinates of stellarator model
def mapCoordinates(u: float, v: float) -> tuple[float, float, float]:
    s = 1.0 # s will always be defined as 1.0
    phi = u * 2 * math.pi # Map u to phi
    theta = v * 2 * math.pi  # Map v to theta
    x, y, z = vmec.vmec2xyz(s, theta, phi) # Now the coordinates can be expressed in cartesian form
    return x, y, z 

# Find the maximum disance in the z direction of the space occupied of the stellarator model.
# In other words, this is the height of the model from the bottom to the top.
# Returns - Ending height (maxz), starting height (minz), and radius (maxr) of stellarator model
def searchMax():
    # Initalize variables for height and radius
    maxz = 0
    minz = 0
    maxr = 0
    s = 1.0 # s will always be defined as 1.0
    # With numpy, the maximumheight, minimum height, and radius can be found
    # Phi and theta are first searched through a period
    # Then the vmec data is mapped to cylindrical coordinates
    # If a new maxz, minz, or maxr is found then it is assigned to the variable(s)
    for phi in np.linspace(0,np.pi/2,num=181):
        for theta in np.linspace(0,np.pi * 2,num=721):
            r, p, z = vmec.vmec2rpz(s, theta, phi)
            if (z > maxz):
                maxz = z
            if (z < minz):
                minz = z
            if (r > maxr):
                maxr = r
    return maxz, minz, maxr

def parametricCurveTest(theta, u):
    s = 1.0  # s will always be defined as 1.0
    phi = u * 2 * math.pi  # Map u to phi
    x, y, z = vmec.vmec2xyz(s, theta, phi)  # Calculate the coordinates
    return x, y, z

def main():
    # The stellarator model can now be made with the .parametricSurface function
    stellarator = cq.Workplane("XY").parametricSurface(mapCoordinates,60,0.0,1.0).val()
    # In order for slicing to occur, the stellarator object needs to be a solid
    # This can be done with OCP

    shell = OCP.TopoDS.TopoDS_Shell()
    builder = OCP.BRep.BRep_Builder()
    builder.MakeShell(shell)
    builder.Add(shell, stellarator.wrapped)
    solid = cq.Solid.makeSolid(cq.Shell(shell))
    storageList = []
    
    # Define parameters for splitting (toroidal)
    slicesz = 10 # Number of slices in the longitudinal direction
    slicesphi = 40 # Number of slices in the latitudinal direction
    maxz, minz, maxr = searchMax()
    thetaIncrement = (2 * math.pi) / slicesphi # The increment for theta (how much the slicing plane will turn) will be 360 degrees divided by the number of slices needed
    totalHeight = maxz - minz  # Total height is the distance between the max and min z values
    heightIncrement = totalHeight / slicesz # The change in height for the slicing plane will be the total height divided by the number of slices needed
    # Initialize cosine and sine values for the longitudinal direction
    tcos = 0.0
    tsin = 0.0
    phi = 0
    theta = 0

    # In the poloidal direction, split stellarator
    for index in range(slicesphi):
        angle = (index * thetaIncrement)
        tcos = math.cos(angle)
        tsin = math.sin(angle)
        splittingPlaneLong = cq.Face.makePlane(None, None, cq.Vector(0, 0, 0), cq.Vector(tcos, tsin, 0))
        stellarator = stellarator.split(splittingPlaneLong)

    # TODO theta 0 to 2pi, splitting with parametricCurve, theta and s constant
    # In the toroidal direction, split stellarator
    for theta in np.linspace(0, 2 * math.pi, num=slicesz, endpoint=False):
        slicingCurve = cq.Workplane("XY").parametricCurve(lambda u: parametricCurveTest(theta, u), 60, 0.0, 1.0).val()
        stellarator = stellarator.split(slicingCurve)
        show_object(slicingCurve)

    storageList = iter(stellarator)

    show_object(stellarator) # Display resulting stellarator model to cq-editor
    export(stellarator) # Call the export function to export model as 3D model file

# Call the main function to start processes
main()
