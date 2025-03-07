
# Extract objects on *all* IHP layers in GDSII file
# Use result with gmsh API directly, no intermediate geo file
# Usage: gds2gmsh <input.gds> 

# File history: 
# Initial version 1 August 2024 Volker Muehlhaus 
# v5: physical groups detection per z position (not per GDSII layer), fragment IHP metal and via layers
# v6: restructured layer data
# v7: add dielectrics
# v8: add physical groups for surfaces 
# v9: add output for *.sif materials, bodies and surfaces
# v10: add conductivity for metal layers and vias
# v11: add surface at substrate group to output
# v12: add mesh size seed, set 0 to disable


import gdspy
import sys
from pathlib import Path
import gmsh

# ============= settings ===========

# enable/disable via array merging
enable_viamerging = True

# oversize of simulation boundary compared to GDSII layout
boundary_distance_xy = 50
boundary_distance_z  = 50

# global factor for via array merging with large spacings
via_oversize_factor = 1

# seed for maximum mesh cell size, set 0 to diable
meshseed = 0

# switch between full IHP SG13G2 stackup and a simple stackup with just one oxide layer (no silicon below, no air above)
fullstack = False

# ============= utilities ==========

def float2string (value):
  return "{:.3f}".format(value)    # fixed 3 decimal digits


# ============= technology specific stuff ===============

# list of purpose to evaluate
purposelist = [
0
]

# Dictionary of names and properties for each GDSII layer number
# Order matters for Vmin and TopVia1, list ist processed top to bottom  
layers = {
 1: { "name": "Activ",      "zmin": 283.7500, "thick": 0.4000, "cond": 35.714E4},
 8: { "name": "Metal1",     "zmin": 284.7900, "thick": 0.4200, "cond": 2.1640E7},
10: { "name": "Metal2",     "zmin": 285.7500, "thick": 0.4900, "cond": 2.3190E7}, 
30: { "name": "Metal3",     "zmin": 286.7800, "thick": 0.4900, "cond": 2.3190E7},
50: { "name": "Metal4",     "zmin": 287.8100, "thick": 0.4900, "cond": 2.3190E7},
67: { "name": "Metal5",     "zmin": 288.8400, "thick": 0.4900, "cond": 2.3190E7},
36: { "name": "MIM",        "zmin": 289.3541, "thick": 0.1500, "cond": 50.000E4},
126:{ "name": "TopMetal1",  "zmin": 290.1803, "thick": 2.0000, "cond": 2.7800E7},
134:{ "name": "TopMetal2",  "zmin": 294.9803, "thick": 3.0000, "cond": 3.0300E7},
6:  { "name": "Cont",       "zmin": 284.1500, "thick": 0.6400, "cond": 2.3900E6, "spacing": 0.5},
19: { "name": "Via1",       "zmin": 285.2100, "thick": 0.5400, "cond": 1.6600E6, "spacing": 0.5},
29: { "name": "Via2",       "zmin": 286.2400, "thick": 0.5400, "cond": 1.6600E6, "spacing": 0.5},
49: { "name": "Via3",       "zmin": 287.2700, "thick": 0.5400, "cond": 1.6600E6, "spacing": 0.5},
66: { "name": "Via4",       "zmin": 288.3000, "thick": 0.5400, "cond": 1.6600E6, "spacing": 0.5},
129:{ "name": "Vmim",       "zmin": 289.5043, "thick": 0.6760, "cond": 2.1910E6, "spacing": 1.5},
125:{ "name": "TopVia1",    "zmin": 289.3300, "thick": 0.8503, "cond": 2.1910E6, "spacing": 1.5},
133:{ "name": "TopVia2",    "zmin": 292.1803, "thick": 2.8000, "cond": 3.1430E6, "spacing": 2.0}
}


# get layername/materialname from GDSII layer number 
def get_metal_layername (num):
  layer = layers.get(num,{"name":"unknown"})
  layername = layer.get("name")
  return layername

# get layer starting position from GDSII layer number 
def get_metal_zmin (num):
  layer = layers.get(num,{"zmin":0})
  zmin = layer.get("zmin")
  return zmin

# get layer thickness position from GDSII layer number 
def get_metal_thickness (num):
  layer = layers.get(num,{"thick":0})
  thickness = layer.get("thick")
  return thickness

# get spacing for via layers, return 0 if not via layer or not in list
def get_via_spacing  (num):
  layer = layers.get(num,{"spacing":0})
  spacing = layer.get("spacing",0)
  return spacing

# get metal conductivity from GDSII layer number 
def get_metal_cond (num):
  layer = layers.get(num,{"cond":0})
  cond = layer.get("cond")
  return cond

# list of dielectric layers

if fullstack:
  dielectrics = [
    {"name":"Substrate",     "zmin": 200.0,  "thick":  80.0,   "eps": 11.9, "sigma": 2},
    {"name":"EPI",           "zmin": 280.0,  "thick":  3.75,  "eps": 11.9, "sigma": 5},
    {"name":"SiO2",          "zmin": 283.75, "thick": 15.73,  "eps":  4.1, "sigma": 0},
  #  {"name":"Passive",       "zmin": 299.48, "thick":  0.4,   "eps":  6.6, "sigma": 0},
    {"name":"Air",           "zmin": 299.48, "thick": 50.0,   "eps":  1.0, "sigma": 0}
  ]
else:
  # simplified list of dielectric layers
  dielectrics = [
    {"name":"SiO2",          "zmin": 270, "thick": 60,  "eps":  4.1, "sigma": 0}
  ]
   


def merge_via_array (polygons, maxspacing):
  # Via array merging consists of 3 steps: oversize, merge, undersize
  # Value for oversize depends on via layer
  # Oversized vias touch if each via is oversized by half spacing
  
  offset = maxspacing/2 + 0.01
  
  offsetpolygonset=gdspy.offset(polygons, offset, join='miter', tolerance=2, precision=0.001, join_first=False, max_points=199)
  mergedpolygonset=gdspy.boolean(offsetpolygonset, None,"or", max_points=199)
  mergedpolygonset=gdspy.offset(mergedpolygonset, -offset, join='miter', tolerance=2, precision=0.001, join_first=False, max_points=199)
  
  # offset and boolean return PolygonSet, we only need the list of polygons from that
  return mergedpolygonset.polygons 


def get_layer_volumes (model, layer):
  # return all volume tags for a given  layer number
 
  # get z position and thickness for this layer
  layer_zpos  = get_metal_zmin(layer)
  layer_thick = get_metal_thickness(layer)

  # get volumes on this layer
  zmin = layer_zpos - 0.01
  zmax = layer_zpos + layer_thick + 0.01
      
  volume_on_layer_list = model.getEntitiesInBoundingBox(-10000,-10000,zmin,10000,10000,zmax,3)
  return  volume_on_layer_list

# ============= main ===============

if len(sys.argv) >= 2:
  input_name = sys.argv[1]
  print ("Input file: ", input_name)
  # get basename of input file, append suffix to identify output polygons
  mesh_name = Path(input_name).stem + ".msh"
   
  # Read GDSII library
  input_library = gdspy.GdsLibrary(infile=input_name)

  # evaluate only first top level cell
  toplevel_cell_list = input_library.top_level()
  cell = toplevel_cell_list[0]
  
  # initialize values for bounding box calculation
  bb_xmin=10000
  bb_ymin=10000
  bb_xmax=-10000
  bb_ymax=-10000

  # list for materials, bodies and surfaces
  sif_materials_list = []
  sif_bodies_list = []
  sif_surfaces_list = []
  
  
  gmsh.initialize()
    
  # iterate over IHP technology layers
  for layer_to_extract in layers.keys():
    
    print ("Evaluating layer ", str(layer_to_extract))
    # flatten hierarchy below this cell
    cell.flatten(single_layer=None, single_datatype=None, single_texttype=None)
    
    # get layers used in cell
    used_layers = cell.get_layers()

    # check if via layer and get max spacing value (used later for via array merging), 
    viaspacing = get_via_spacing(layer_to_extract)  # value is 0 if not via layer
    if (viaspacing>0) and enable_viamerging:
      print ("Detected via layer, max via spacing in array = ", str(viaspacing))
    
        
    # check if layer-to-extract is used in cell 
    if (layer_to_extract in used_layers):
            
      # iterate over layer-purpose pairs (by_spec=true)
      # do not descend into cell references (depth=0)
      LPPpolylist = cell.get_polygons(by_spec=True, depth=0)
      for LPP in LPPpolylist:
        layer = LPP[0]
        purpose = LPP[1]
        
        # now get polygons for this one layer-purpose-pair
        if (layer==layer_to_extract) and (purpose in purposelist):
          
          layername = get_metal_layername(layer)
          layerpolygons = LPPpolylist[(layer, purpose)]

          # for via layers do via array merging
          if (viaspacing>0) and enable_viamerging:
            layerpolygons = merge_via_array (layerpolygons, viaspacing*via_oversize_factor)
                    
          # get z position and thickness for this layer
          layer_zpos  =get_metal_zmin(layer)
          layer_thick =get_metal_thickness(layer)
                    
          # iterate over layer polygons
          for polypoints in layerpolygons:
            
            # write polygon info to comment line
            numvertices = int(polypoints.size/polypoints.ndim)

            # empty list of boundary line tags
            linetaglist = []
            vertextaglist = []
            
            # define vertices
            for vertex in range(numvertices):
              
              x = polypoints[vertex,0]
              y = polypoints[vertex,1]
              
              # addPoint parameters: x (double), y (double), z (double), meshSize = 0. (double), tag = -1 (integer)
              meshsize = meshseed
              vertextag = gmsh.model.occ.addPoint(x,y,layer_zpos,meshsize,-1)
              vertextaglist.append(vertextag)
                            
              # update bounding box information
              if x<bb_xmin: bb_xmin=x
              if x>bb_xmax: bb_xmax=x
              if y<bb_ymin: bb_ymin=y
              if y>bb_ymax: bb_ymax=y
           
            # after writing the vertices, we combine them to boundary lines
            for v in range(numvertices):
              pt_start = vertextaglist[v]
              if v==(numvertices-1):
                pt_end = vertextaglist[0]
              else:
                pt_end = vertextaglist[v+1]

              # addLine parameters: startTag (integer), endTag (integer), tag = -1 (integer)
              linetag = gmsh.model.occ.addLine(pt_start, pt_end, -1)
              linetaglist.append(linetag)
              

            # after creating the lines, we can create a curve loop and a surface 
            # for that, we need the line segment numbers again
            curvetag   = gmsh.model.occ.addCurveLoop(linetaglist, tag=-1)
            surfacetag = gmsh.model.occ.addPlaneSurface([curvetag], tag=-1)
            # print('Created curve   with tag ' + str(curvetag))
            # print('Created surface with tag ' + str(surfacetag))
            
            returnval = gmsh.model.occ.extrude([(2,surfacetag)],0,0,layer_thick)
            # print('Return value frome extrude ' + str(returnval))
            volumetag = returnval[1][1]
            # print('Volume tag ' + str(volumetag))

            gmsh.model.occ.synchronize()



  # We have created initial 3D volumes from GDSII, now iterate over 3D entities to merge them

  volumelist = gmsh.model.getEntities(3)
  volumecount = len(volumelist)
  if volumecount>0:
    print('Number of volumes = ' + str(volumecount)) 

    # iterate over IHP technology layers again and postprocess layer objects

    # try to merge volumes on each layer
    print ("\nBegin Merging")


    for layer in layers.keys():
      volume_on_layer_list = get_layer_volumes(gmsh.model, layer)

      # try boolean union of volumes on this layer
      if len(volume_on_layer_list)>1:
        # get first element and delete from list
        first = volume_on_layer_list.pop(0)
        print('  Layer = ' + get_metal_layername(layer)) 
        print('  FUSE, object = ' + str(first)) 
        print('  FUSE, tool   = ' + str(volume_on_layer_list)) 
        
        gmsh.model.occ.fuse([first],volume_on_layer_list, -1)
        gmsh.model.occ.synchronize()

    
    # fragment adjacent layers with the vias, to create shared surfaces which are necessary for meshing

    # fragment Activ, Cont and Metal1
    Activ_list   = get_layer_volumes(gmsh.model,1)
    Cont_list    = get_layer_volumes(gmsh.model,6)
    Metal1_list  = get_layer_volumes(gmsh.model,8)

    gmsh.model.occ.fragment(Activ_list,Cont_list)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.fragment(Cont_list,Metal1_list)
    gmsh.model.occ.synchronize()

    Via1_list    = get_layer_volumes(gmsh.model,19)
    Metal2_list  = get_layer_volumes(gmsh.model,10)

    gmsh.model.occ.fragment(Metal1_list,Via1_list)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.fragment(Via1_list,Metal2_list)
    gmsh.model.occ.synchronize()

    Via2_list    = get_layer_volumes(gmsh.model,29)
    Metal3_list  = get_layer_volumes(gmsh.model,30)

    gmsh.model.occ.fragment(Metal2_list,Via2_list)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.fragment(Via2_list,Metal3_list)
    gmsh.model.occ.synchronize()

    Via3_list    = get_layer_volumes(gmsh.model,49)
    Metal4_list  = get_layer_volumes(gmsh.model,50)

    gmsh.model.occ.fragment(Metal3_list,Via3_list)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.fragment(Via3_list,Metal4_list)
    gmsh.model.occ.synchronize()

    Via4_list    = get_layer_volumes(gmsh.model,66)
    Metal5_list  = get_layer_volumes(gmsh.model,67)

    gmsh.model.occ.fragment(Metal4_list,Via4_list)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.fragment(Via4_list,Metal5_list)
    gmsh.model.occ.synchronize()

    TopVia1_list     = get_layer_volumes(gmsh.model,125)
    TopMetal1_list   = get_layer_volumes(gmsh.model,126)

    gmsh.model.occ.fragment(Metal5_list,TopVia1_list)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.fragment(TopVia1_list,TopMetal1_list)
    gmsh.model.occ.synchronize()

    Vmim_list       = get_layer_volumes(gmsh.model,129)
    MIM_list        = get_layer_volumes(gmsh.model,36)
    
    gmsh.model.occ.fragment(MIM_list,Vmim_list)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.fragment(Vmim_list,TopMetal1_list)
    gmsh.model.occ.synchronize()

    TopVia2_list   = get_layer_volumes(gmsh.model,133)
    TopMetal2_list = get_layer_volumes(gmsh.model,134)

    gmsh.model.occ.fragment(TopMetal1_list,TopVia2_list)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.fragment(TopMetal2_list,TopVia2_list)
    gmsh.model.occ.synchronize()
        

    
    # assign physical groups for metals
    # layers are identified by z position    
    # special case Vmim: this overlaps in z-direction with TopVia1, so we need to process Vmim first
    Vmim_volume_tags = []

    print('Setting physical groups and names for metals')
    # iterate over layer numbers
    for layer in layers.keys():
      volume_on_layer_list = get_layer_volumes(gmsh.model,layer)
      layername = get_metal_layername(layer)
      layer_volume_tags = []  
      
      volume_index = 0 # counter needed for naming the surfaces
      for volume in volume_on_layer_list:
        volume_index = volume_index+1

        if layername=='Vmim':
          Vmim_volume_tags.append (volume[1])

        if ((layername=='TopVia1') and (volume[1] in Vmim_volume_tags)):
          print('  volume = ' + str(volume[1]) + ' skip, already registered for layer Vmim') 
        else:
          # we have a valid volume on this layer
          volumetag = volume[1]
          layer_volume_tags.append (volumetag)
          print('  volume = ' + str(volumetag) + ' found on layer ' + layername) 

          # We have found a volume, now get the surfaces of this volume, so that we can assign a physical group for the surfaces also
          surfacelooptags, surfacetags = gmsh.model.occ.getSurfaceLoops(volumetag)
          for surfacetag in surfacetags:
            surfacename = layername + '_' + str(volume_index)
            physical_tag = gmsh.model.addPhysicalGroup(2, surfacetag, tag=-1)
            gmsh.model.setPhysicalName(2, physical_tag, surfacename) 
            #print('    surface = ' + str(surfacename) + '  ' + str(surfacetag))
          
      
      physical_tag = gmsh.model.addPhysicalGroup(3, layer_volume_tags, tag=-1)
      gmsh.model.setPhysicalName(3, physical_tag, layername) 
      
      # add material and body to lists
      # for now, we define metals as epsr=1
      if volume_index > 0:
        cond = get_metal_cond(layer)
        material_index = len(sif_materials_list) + 1  # index of item that we add now
        sif_materials_list.append ('Material ' + str(material_index) +'\n  Name = "' + layername + '"\n  Relative Permittivity = 1\n  Relative Permeability = 1.0\n  Electric Conductivity = ' + str(cond) + '\nEnd\n\n' )

        body_index = len(sif_bodies_list) + 1  # index of item that we add now
        sif_bodies_list.append ('Body ' + str(body_index) +'\n  Name = "' + layername + '"\n  Material = ' + str(material_index) +'\n  Equation = 1\nEnd\n\n' )

            


    # calculate simulation boundary from bounding box of GDSII polygons

    boundary_xmin = bb_xmin - boundary_distance_xy
    boundary_xmax = bb_xmax + boundary_distance_xy
    boundary_ymin = bb_ymin - boundary_distance_xy
    boundary_ymax = bb_ymax + boundary_distance_xy
    dx = boundary_xmax - boundary_xmin
    dy = boundary_ymax - boundary_ymin
    
    # add substrate and dielectric layers
    dieelectric_volume_tags = []

    for dielectric in dielectrics:
      zmin = dielectric.get("zmin", 0) 
      dz   = dielectric.get("thick",1) 
      name = dielectric.get("name", "unnamed")

      # increase thickness of air layer
      if name=="air":
        dz = dz + boundary_distance_z

      # dielectric box
      box = gmsh.model.occ.addBox(boundary_xmin, boundary_ymin,zmin, dx, dy, dz, -1)
      dieelectric_volume_tags.append(box)

    gmsh.model.occ.synchronize()

    # fragment dieelectrics to each other to align the meshes
    print('Fragment dielectrics with dielectrics')
    for i in range(len(dieelectric_volume_tags)-1):
      tag = dieelectric_volume_tags[i]   
      a,b = gmsh.model.occ.fragment([(3, tag)],[(3,tag+1)])
            
    gmsh.model.occ.synchronize()

    # fragment metals to dielectrics inside their volume, cut out metals from dielectrics
    print('Fragment dielectrics with metals')
    for i in range(len(dieelectric_volume_tags)):
       dielectric_tag = dieelectric_volume_tags[i] 

       # get bounding box of dielectric
       dk_xmin, dk_ymin, dk_zmin, dk_xmax, dk_ymax, dk_zmax = gmsh.model.getBoundingBox(3, dielectric_tag)
       # get metals inside that box, but this also includes the dielectrit itself
       included_volumes = gmsh.model.getEntitiesInBoundingBox(dk_xmin, dk_ymin,dk_zmin-0.01,dk_xmax,dk_ymax,dk_zmax+0.01, 3)  
       if len(included_volumes)>1:
         a,b = gmsh.model.occ.fragment([(3, dielectric_tag)],included_volumes)

    gmsh.model.occ.synchronize()

    # now dielectric are fully fragmented, identify them again by comparing to the z position from the list
    print('Setting physical groups and names for dielectrics')
    for dielectric in dielectrics:
      zmin = dielectric.get("zmin", 0) 
      zmax = zmin +  dielectric.get("thick",1) 
      name = dielectric.get("name", "unnamed")
      eps  = dielectric.get("eps",1) 
      cond = dielectric.get("cond",0) 

      
      # get all volumes in that layer
      all_volumes = gmsh.model.getEntitiesInBoundingBox(-10000,-10000,zmin-0.1,10000,10000,zmax+0.1,3)
      # get metals in that layer, so that we can exclude them
      metal_volumes = gmsh.model.getEntitiesInBoundingBox(bb_xmin,bb_ymin,zmin-0.1,bb_xmax,bb_ymax,zmax+0.1,3)
      for metal in metal_volumes:
        all_volumes.remove(metal)
      
      dk_volume_tags = [model[1] for model in all_volumes]
      #print('  dk_volume_tags = ' + str(dk_volume_tags))
      physical_tag = gmsh.model.addPhysicalGroup(3, dk_volume_tags, tag=-1)
      #print('  physical_tag = ' + str(physical_tag))
      gmsh.model.setPhysicalName(3, physical_tag, name) 

      # add material and body to lists
      material_index = len(sif_materials_list) + 1  # index of item that we add now
      sif_materials_list.append ('Material ' + str(material_index) +'\n  Name = "' + name + '"\n  Relative Permittivity = ' + str(eps) +'\n  Relative Permeability = 1\n  Electric Conductivity = ' + str(cond) + '\nEnd\n\n' )

      body_index = len(sif_bodies_list) + 1  # index of item that we add now
      sif_bodies_list.append ('Body ' + str(body_index) +'\n  Name = "' + name + '"\n  Material = ' + str(material_index) +'\n  Equation = 1\nEnd\n\n' )


    # create physical group for surfaces at outer simulation boundary
    # get model min and max z values
    boundary_zmin =  1000
    boundary_zmax = -1000
    for dielectric in dielectrics:
      zmin = dielectric.get("zmin", 0) 
      zmax = zmin +  dielectric.get("thick",1) 
      if zmin < boundary_zmin:
        boundary_zmin=zmin
      if zmax > boundary_zmax:
        boundary_zmax = zmax  

    # get all surfaces at xmin
    delta = 0.1
    surfaces_at_xmin = gmsh.model.getEntitiesInBoundingBox(boundary_xmin-delta, boundary_ymin-delta,boundary_zmin-delta, boundary_xmin+delta, boundary_ymax+delta,boundary_zmax+delta, 2 )
    surfaces_at_xmax = gmsh.model.getEntitiesInBoundingBox(boundary_xmax-delta, boundary_ymin-delta,boundary_zmin-delta, boundary_xmax+delta, boundary_ymax+delta,boundary_zmax+delta, 2 )
    surfaces_at_ymin = gmsh.model.getEntitiesInBoundingBox(boundary_xmin-delta, boundary_ymin-delta,boundary_zmin-delta, boundary_xmax+delta, boundary_ymin+delta,boundary_zmax+delta, 2 )
    surfaces_at_ymax = gmsh.model.getEntitiesInBoundingBox(boundary_xmin-delta, boundary_ymax-delta,boundary_zmin-delta, boundary_xmax+delta, boundary_ymax+delta,boundary_zmax+delta, 2 )
    surfaces_at_zmin = gmsh.model.getEntitiesInBoundingBox(boundary_xmin-delta, boundary_ymin-delta,boundary_zmin-delta, boundary_xmax+delta, boundary_ymax+delta,boundary_zmin+delta, 2 )
    surfaces_at_zmax = gmsh.model.getEntitiesInBoundingBox(boundary_xmin-delta, boundary_ymin-delta,boundary_zmax-delta, boundary_xmax+delta, boundary_ymax+delta,boundary_zmax+delta, 2 )
    
    outer_surfaces = surfaces_at_xmin + surfaces_at_xmax + surfaces_at_ymin + surfaces_at_ymax + surfaces_at_zmin + surfaces_at_zmax
    # print('  outer_surfaces = ' + str(outer_surfaces))

    outer_surface_tags = [model[1] for model in outer_surfaces]
    physical_tag = gmsh.model.addPhysicalGroup(2, outer_surface_tags, tag=-1)
    gmsh.model.setPhysicalName(2, physical_tag, 'Model_boundary') 

    # create physical group for surfaces at substrate (top of EPI)
    if len(dielectrics)>1:
      dielectric = dielectrics[1]
      assert dielectric["name"]=="EPI"
      EPI_zmax = dielectric.get("zmin", 0) + dielectric.get("thick",1) 
      surfaces_at_EPI = gmsh.model.getEntitiesInBoundingBox(boundary_xmin-delta, boundary_ymin-delta,EPI_zmax-delta, boundary_xmax+delta, boundary_ymax+delta,EPI_zmax+delta, 2 )

      print('  surfaces_at_EPI = ' + str(surfaces_at_EPI))

      EPI_surface_tags = [model[1] for model in surfaces_at_EPI]
      physical_tag = gmsh.model.addPhysicalGroup(2, EPI_surface_tags, tag=-1)
      gmsh.model.setPhysicalName(2, physical_tag, 'EPI_surface') 



  gmsh.model.occ.synchronize()
 
  # Launch the GUI to see the model:
  gmsh.fltk.run()

    # write plain geo with no mesh
  geo_name = Path(input_name).stem + ".geo_unrolled"
  gmsh.write(geo_name)


  # generate mesh and show again
  gmsh.option.setNumber("Mesh.SaveAll", 1)  # save everything, no matter if in physical group or not
  gmsh.model.mesh.generate()
  gmsh.fltk.run()
  gmsh.write(mesh_name)


  # We can use this to clear all the model data:
  gmsh.clear()
  gmsh.finalize()

  # output file for Elmer *.sif structure information file
  sif_name = Path(input_name).stem + ".sif"
  sif_file = open(sif_name, "w")

  for item in sif_materials_list:
    sif_file.write(item)
  
  for item in sif_bodies_list:
    sif_file.write(item)

  for item in sif_surfaces_list:
    sif_file.write(item)


  sif_file.close()

  
else:
  print ("Usage: gds2gmsh <input.gds> ")

  
  
