## Forest inventory processing
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script loads the masked SDMs prepared in script "3-sdm-mask", creates 
## binary maps by thresholding probability values >= 0.6, and sums all species
## to obtain a map of SDM overlap expressed in number of species. Lastly, it 
## calculates area covered by species' distributions. 
##
## Outputs:
## - tif, binary SDM maps for each selected species
## - tif, SDM overlap across all selected species

import fiona
from shapely.geometry import box, mapping
from fiona.crs import from_epsg
import rasterio

# Function to split bounding box into tiles
def split_bbox(bbox, rows, cols):
    min_x, min_y, max_x, max_y = bbox
    width = (max_x - min_x) / cols
    height = (max_y - min_y) / rows
    
    tiles = []
    for i in range(rows):
        for j in range(cols):
            tile_min_x = min_x + (j * width)
            tile_min_y = min_y + (i * height)
            tile_max_x = tile_min_x + width
            tile_max_y = tile_min_y + height
            tiles.append((tile_min_x, tile_min_y, tile_max_x, tile_max_y))
    return tiles

# Set bounding box
model_bbox = (-942562.735316926, 5971.82034726957, 429467.232992499, 744578.222380984)

# Get tiles
tiles = split_bbox(model_bbox, 10, 8)

# Stream-filter the vector layer by bounding box
for i in range(len(tiles)):
    with fiona.open("Data/Spatial/Ecoforestry/CARTE_ECO_MAJ_PROV.gpkg", layer="pee_maj") as src:
        bbox = tiles[i]
        
        schema = src.schema
        crs = src.crs
        
        with fiona.open((f"Intermediate/Ecoforestry/ecoforestry_tile_{i}.gpkg"), "w", driver="GPKG", schema=schema, crs=crs) as dst:
            for feature in src.filter(bbox=bbox):
                dst.write(feature)
