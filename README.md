# gds2gmsh

This tool is used in the IHP SG13G2 RFIC design flow, to convert layout
in GDSII file format to a model/mesh for the gmsh meshing tool.
https://gmsh.info/
In addition to creating the gmsh files, material properties are stored 
in *.sif format for use with the ElmerFEM solver.

# Prerequisites

To run the tool, Python libraries gdspy and gmsh must be installed.

# Usage

To build the gmsh model/mesh from the GDSII input file, specify the 
*.gds input filename as the first parameter.

example:
```
python gds2gmsh.py ind_square_100x100.gds
```

This will create 3 output files:
the *.msh file with the mesh from gmsh (no gemetry), 
the *.geo_unrolled file with the geometry for gmsh (no mesh) and
an additional *.sif file with material information for use with ElmerFEM.

The gmsh user interface is started to show the geometry. 

![plot](./doc/geometry_view.png)

When this window is closed, the user interface comes up again, 
showing the resulting mesh.

Mesh creation is controlled by parameter 'meshseed' that is defined in the code, 
and might need further improvement.

![plot](./doc/mesh_view.png)




# Theory of operation

The core of layout conversion from GDSII to gmsh is to flatten the layout 
into a flat list of polygons, iterate over all vertices of each polygon and 
build a list of points (1-D) which are then combined to lines of adjacent points (2-D)
which are then combined to a curve loop (2-D) which is then declared as surface (2-D).
For each 2-D surface, the z-position is determined from the source layer by a lookup 
table with stackup information, i.e. z positions of each layer. 
The 2-D surfaces are then extruded in z-direction to give them a finite height, 
as defined in the stackup lookup table, and a 3-D volume object is created from 
the extruded polygons.

For the surrounding dielectric layers (oxide, passivation) and the semiconductor 
substrate below, there is no layout data in the GDSII file, so the xy size of 
these dielectrics are determined from the bounding box of the layer’s layout data, 
plus some margin. Again, the z-position and thickness are obtained from the stackup
lookup table.

Physical 3D groups are defined for each metal layer or dielectric, 
and physical 2D surfaces are defined for each metal or via polygon.









