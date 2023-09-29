# Syn Hubbard (CNERG)

# Import statements for calculations and export
import cadquery as cq
import math
from cadquery import exporters
import read_vmec
import OCP
import numpy as np
vmec = read_vmec.vmec_data(r'C:\Users\vich2\Downloads\plas_eq.nc') # Global scope of data

# Export function
def export(surface):
    exporters.export(surface, './stellarator.step')

# Using the data, lets map the coordinates
def func(u: float, v: float) -> tuple[float, float, float]:
    s = 1.0
    phi = u * 2 * math.pi # Map u to phi
    theta = v * 2 * math.pi  # Map v to theta
    x, y, z = vmec.vmec2xyz(s, theta, phi)
    return x, y, z

def searchMax():
    maxz = 0
    minz = 0
    maxr = 0
    s = 1
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

def main():
    stellarator = cq.Workplane("XY").parametricSurface(func,60,0.0,1.0).val()
    shell = OCP.TopoDS.TopoDS_Shell()
    builder = OCP.BRep.BRep_Builder()
    builder.MakeShell(shell)
    builder.Add(shell, stellarator.wrapped)
    solid = cq.Solid.makeSolid(cq.Shell(shell))
    
    # Define parameters for splitting, in reality, it would be nice for this to be inputted by the user
    slicesz = 10  # replace integer literals -> prompt for slices
    slicesphi = 40
    maxz, minz, maxr = searchMax()
    shellThickness = 0.01
    thetaIncrement = (2 * math.pi) / slicesphi
    totalHeight = maxz - minz  # Total height of the sphere
    heightIncrement = totalHeight / slicesz
    Tcos = 0.0
    Tsin = 0.0

    # In the longitudinal direction, split sphere
    for index in range(slicesphi):
        angle = (index * thetaIncrement)
        Tcos = math.cos(angle)
        Tsin = math.sin(angle)
        splittingPlaneLong = cq.Face.makePlane(None, None, cq.Vector(0, 0, 0), cq.Vector(Tcos, Tsin, 0))
        stellarator = stellarator.split(splittingPlaneLong)

    # In the latitudinal direction, split sphere with evenly spaced slices
    for index in range(slicesz):  # Start from 1 to avoid splitting at the top and bottom poles
        height = minz + (index * heightIncrement)  # Adjust height based on the increment
        if (height == 0):
            continue
        splittingPlaneLat = cq.Face.makePlane(None, None, cq.Vector(0, 0, height), cq.Vector(0, 0, 1))
        stellarator = stellarator.split(splittingPlaneLat)

    show_object(stellarator)
    export(stellarator)

main()
