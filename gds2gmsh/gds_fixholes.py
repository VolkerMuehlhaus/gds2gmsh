# Extract objects on *all* IHP layers in GDSII file
# Fix  polygons with holes by slicing them into separate polygons
# Usage: gds_fixholes <input.gds> 

# File history: 
# Initial version 3 August 2024 Volker Muehlhaus 

import gdspy
import sys
from pathlib import Path

# ============= technology specific stuff ===============

# list of layers to evaluate, we only check  metal layers

layerlist = [
1,
8,
10,
30,
50,
67,
126,
134
]

# list of purpose to evaluate
purposelist = [
0,
70
]


# ============= utilities ==========

def float2string (value):
  return "{:.3f}".format(value)    # fixed 3 decimal digits



# ============= main ===============

if len(sys.argv) >= 2:
  input_name = sys.argv[1]
  print ("Input file: ", input_name)
  # get basename of input file, append suffix to identify output polygons
  output_name = Path(input_name).stem + "_sliced.gds"
    
  # Read GDSII library
  output_library = gdspy.GdsLibrary(infile=input_name)
  
  # iterate over cells
  for cell in output_library:
    print('\ncellname = ' + str(cell.name))

    # iterate over polygons
    for poly in cell.polygons:
      
      if len(poly.polygons)>0:
        # points of this polygon
        polypoints = poly.polygons[0]

        poly_layer = poly.layers[0]
        poly_purpose = poly.datatypes[0]

        if ((poly_layer in layerlist) and (poly_purpose in purposelist)):
        
          # get number of vertices
          numvertices = len(polypoints) 
          
          seen   = set()    # already seen vertex values
          dupefound = False

          # iterate over vertices to find duplicates
          for i_vertex in range(numvertices):
            
            # print('polypoints  = ' + str(polypoints))
            x = polypoints[i_vertex][0]
            y = polypoints[i_vertex][1]
            
            # create string representation so that we can check for duplicates
            vertex_string = str(x)+','+str(y)
            if vertex_string in seen:
              dupefound = True
              print('      found duplicate at vertex ' + str(i_vertex) + ': ' + vertex_string)
            else:
              seen.add(vertex_string)  

          if dupefound:
                        
            # do the slicing
            
            # convert polygon to format required for slicing
            basepoly_points = []

            for i_vertex in range(numvertices):
              basepoly_points.append((polypoints[i_vertex,0], polypoints[i_vertex,1]))

            # create new polygon
            basepoly = gdspy.Polygon(basepoly_points, layer=poly_layer, datatype=poly_purpose)  
            fractured = basepoly.fracture(max_points=6)

            # add fractured polygon to cell
            cell.add(fractured)

            # invalidate original polygon
            poly.layers=[0]
            # remove original polygon
            cell.remove_polygons(lambda pts, layer, datatype:
              layer == 0)
            

  
  
  # write to output file
  output_library.write_gds(output_name)

  
else:
  print ("Usage: gds_fixholes <input.gds> ")

  
  
