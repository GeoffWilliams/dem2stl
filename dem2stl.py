#!/usr/bin/env python
#    
# DEM (Digital Elevation Model) to STL (3d printer)
# Copyright (C) 2014 Geoff Williams<geoff@geoffwilliams.me.uk>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import argparse
import Image
import numpy
import logging
import math


R_RATIO=0.002
G_RATIO=0.001
B_RATIO=-0.002
PIXEL_RATIO=1
RESAMPLE_SIZE=1
RESAMPLE_STEP=0
RESAMPLE_SQUARED=0
MODEL_DEPTH=6
Z_MIN=0
Z_MAX=0

def main():
    global R_RATIO
    global G_RATIO
    global B_RATIO
    global PIXEL_RATIO
    global RESAMPLE_SIZE
    global RESAMPLE_STEP
    global RESAMPLE_SQUARED
    global MODEL_DEPTH

    logging.info("DEM TO STL")
    parser = argparse.ArgumentParser(description='DEM (Digital Elevation Model) to STL (3D printer) conversion')
    # input file
    parser.add_argument('--input-file', dest='input_file', action='store',
            help='filename for input file', required=True)
    # output file
    parser.add_argument('--output-file', dest='output_file', action='store',
            help='filename for output file', required=True)
    
    # r ratio
    parser.add_argument('--r-ratio', dest='r_ratio', action='store',
            help='ratio for red pixels', default=R_RATIO, type=float)
    
    # g ratio
    parser.add_argument('--g-ratio', dest='g_ratio', action='store',
            help='ratio for green pixels', default=G_RATIO, type=float)
    
    # b ratio
    parser.add_argument('--b-ratio', dest='b_ratio', action='store',
            help='ratio for blue pixels', default=B_RATIO, type=float)
    
    # pixel ratio
    parser.add_argument('--pixel-ratio', dest='pixel_ratio', action='store',
            help='pixel to world unit ratio', default=PIXEL_RATIO, type=float)
    
    # resample size
    parser.add_argument('--resample-size', dest='resample_size', action='store',
            help='resample the source image by averaging squares of this size (px)', type=int, default=RESAMPLE_SIZE)
    
    # model depth
    parser.add_argument('--model-depth', dest='model_depth', action='store',
            help='minimum depth of model', type=int, default=MODEL_DEPTH)

    # debug mode
    parser.add_argument('--debug', dest='debug', action='store_true',
            help='enable debug messages', default=False)

    args = parser.parse_args()
    R_RATIO = args.r_ratio
    G_RATIO = args.g_ratio
    B_RATIO = args.b_ratio
    PIXEL_RATIO = args.pixel_ratio
    RESAMPLE_SIZE = args.resample_size    
    RESAMPLE_STEP = RESAMPLE_SIZE * 2 + 1
    RESAMPLE_SQUARED = math.pow(RESAMPLE_STEP,2)
    MODEL_DEPTH = args.model_depth

    if (args.debug):
        print("enabling debug mode...")
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)
        logging.debug("...debug mode enabled")

    logging.debug("RGB and pixel scale ratios:  r=%f g=%f b=%f p=%f" % 
        (R_RATIO, G_RATIO, B_RATIO, PIXEL_RATIO)) 

    convert_height_map(
        height_map(args.input_file), 
        args.output_file)
       
 
       

def height_map(input_file):
    """
    Convert the image pixel values to a 2D array of heights



    pixel order JPEG:

    0,0     --------->X+
    |
    |
    | y+
    v
    """
    global Z_MIN
    global Z_MAX

    im = Image.open(input_file) #Can be many different formats.
    pix = im.load()

    x_dim, y_dim = im.size
    logging.debug("creating map %dx%d" % (x_dim, y_dim))
    values = numpy.empty((x_dim, y_dim))

    # Store the min and max heights as we scan through the image.  This way
    # we can put the output into a base so that negative heights don't go through
    # the base, etc
 
    for x in range(0,x_dim):
        for y in range(y_dim -1,0,-1):
            (r,g,b) = pix[x,y]
            
            # r+g+b = meters * scale factor(???)
            # flip y axis
            z = (r * R_RATIO)+(g * G_RATIO)+(b * B_RATIO)
            Z_MIN = min(Z_MIN, z)
            Z_MAX = max(Z_MAX, z)
            values[x][y_dim-y] = z

    if (RESAMPLE_SIZE > 0):
        values = resample(values)

    logging.debug("min height %f max height %f" % (Z_MIN, Z_MAX))

    return values

def average(values, source_x_dim, source_y_dim, resampled_x, resampled_y):
    """
     resample_size indicates the border around each pixel that will be averaged
     eg for a resample_size of 2, the illustration below shows the pixels that 
     will be averaged for pixel P

     X X X X X
     X X X X X
     X X P X X
     X X X X X
     X X X X X
    """
    
    source_x = resampled_x * RESAMPLE_STEP
    source_y = resampled_y * RESAMPLE_STEP

    # find the x,y values for each corner of the square illustrated above
    x_min = source_x
    x_max = source_x + RESAMPLE_STEP
    y_min = source_y
    y_max = source_y + RESAMPLE_STEP
    
    # take the mean of the pixels indicated by resample size
    v=0
    for x in range(x_min, x_max):
        for y in range(y_min, y_max):
            v += values[x][y]     

    v = v / RESAMPLE_SQUARED
    return v


def resample(values):
    source_x_dim = len(values)
    source_y_dim = len(values[0])

    x_dim = int(round(source_x_dim / RESAMPLE_STEP))
    y_dim = int(round(source_y_dim / RESAMPLE_STEP))

    logging.debug("resample to %f x %f" % (x_dim, y_dim))

    resampled = numpy.empty((x_dim, y_dim))
    for x in range(0, x_dim):
        for y in range(0, y_dim):
            resampled[x][y]=average(values, source_x_dim, source_y_dim, x, y)

    return resampled

def stl_header(output_file):
    output_file.write("solid DEM\n")
    
def stl_footer(output_file):
    output_file.write("endsolid\n")


def convert_height_map(height_map, output_filename):
    """
    Convert height map to STL data


    STL axis order:


    ^ 
    |Y+
    |
    0,0,0 -------> x+


    ( y axis flipped vs jpeg)
    """
    output_file = open(output_filename, "w")

    stl_header(output_file)

    # process each inner pixel leaving a 1px border of unprocessed pixels
    max_x = len(height_map) -1
    max_y = len(height_map[0]) - 1

    #max_x = 3
    #max_y = 3

    for x in range(1, max_x):
        for y in range(1, max_y):
            link_pixel(output_file, height_map, x, y)

    # link the N, S, E, W sides to zero height (the previously unprocessed border)
    stitch_base(output_file, height_map)
    # N
    

    # S

    # E

    # W

    # link the bottom sides to create a solid
    

    stl_footer(output_file)

def stitch_base(output_file, height_map):
    """
    Stich the base by drawing triangles anchoring the top to the base
    using a double triangle for each n - 1 pair of vertices
    .  

       .  

          .
    
    .  .  .

    """

    base = 0 - MODEL_DEPTH - Z_MIN   

    x_min = 0    
    x_max = len(height_map) -1 
    y_min = 0
    y_max = len(height_map[0]) -1


    sx_max = scale_pixel(x_max)
    sy_max = scale_pixel(y_max)

    # N and S (4x triangles)
    for x in range (1, len(height_map)):
        sxm1 = scale_pixel(x-1)
        sx = scale_pixel(x)
 
        # N
        facet(output_file, {
            "xs": [sx, sx, sxm1],
            "ys": [sy_max, sy_max, sy_max],
            "zs": [height_map[x][y_max], base, base],
        })
        facet(output_file, {
            "xs": [sx, sxm1, sxm1],
            "ys": [sy_max, sy_max, sy_max],
            "zs": [height_map[x][y_max], base, height_map[x-1][y_max]],
        })

        # S
        facet(output_file, {
            "xs": [sx, sxm1, sxm1],
            "ys": [y_min, y_min, y_min],
            "zs": [height_map[x][y_min], height_map[x-1][y_min], base]
        })
        facet(output_file, {
            "xs": [sx, sxm1, sx],
            "ys": [y_min, y_min, y_min],
            "zs": [height_map[x][y_min], base, base],
        })


    
    # W and E
    for y in range (1, len(height_map[0])):

        sym1 = scale_pixel(y-1)
        sy = scale_pixel(y)

        # W
        facet(output_file, {
            "xs": [x_min, x_min, x_min],
            "ys": [sy, sy, sym1],
            "zs": [height_map[x_min][y], base, base],
        })

        facet(output_file, {
            "xs": [x_min, x_min, x_min],
            "ys": [sy, sym1, sym1],
            "zs": [height_map[x_min][y], base, height_map[x_min][y-1]],
        })

        # E
        facet(output_file, {
            "xs": [sx_max, sx_max, sx_max],
            "ys": [sy, sym1, sym1],
            "zs": [height_map[x_max][y], height_map[x_max][y-1], base],
        })
        facet(output_file, {
            "xs": [sx_max, sx_max, sx_max],
            "ys": [sy, sym1, sy],
            "zs": [height_map[x_max][y], base, base],
        })

    # close the shape with two triangles on the bottom
    facet(output_file, {
        "xs": [x_min, x_min, sx_max],
        "ys": [y_min, sy_max, sy_max],
        "zs": [base, base, base],
    })

    facet(output_file, {
        "xs": [x_min, sx_max, sx_max],
        "ys": [y_min, sy_max, y_min],
        "zs": [base, base, base],
    })



def facet(output_file, triangle):
    """
    write a facet using arrays of coordinates
    """

     
    # calculate triangle normals - see:  http://www.mathsisfun.com/algebra/vectors-cross-product.html
    ax = triangle["xs"][0] - triangle["xs"][2]
    ay = triangle["ys"][0] - triangle["xs"][2]
    az = triangle["zs"][0] - triangle["xs"][2]

    bx = triangle["xs"][0] - triangle["xs"][1]
    by = triangle["ys"][0] - triangle["ys"][1]
    bz = triangle["zs"][0] - triangle["zs"][1]

    nx = (ay * bz) - (az * by)
    ny = (az * bx) - (ax * bz)
    nz = (ax * by) - (ay * bx)

    # extra credit - normal should be of unit length..


    output_file.write("facet normal %f %f %f\n" % (nx, ny, nz))
    output_file.write("\touter loop\n")
    output_file.write("\t\tvertex %f %f %f\n" % (triangle["xs"][0], triangle["ys"][0], triangle["zs"][0]))
    output_file.write("\t\tvertex %f %f %f\n" % (triangle["xs"][1], triangle["ys"][1], triangle["zs"][1]))
    output_file.write("\t\tvertex %f %f %f\n" % (triangle["xs"][2], triangle["ys"][2], triangle["zs"][2]))
    output_file.write("\tendloop\n")
    output_file.write("endfacet\n")


def triangle(coords, k1, k2, k3):
    """
    construct a map of ordered coordinates
    """
    return {
        "xs": [coords[k1]["x"], coords[k2]["x"], coords[k3]["x"]], 
        "ys": [coords[k1]["y"], coords[k2]["y"], coords[k3]["y"]],
        "zs": [coords[k1]["z"], coords[k2]["z"], coords[k3]["z"]],
    }

def scale_pixel(n):
    return n * PIXEL_RATIO

def link_pixel(output_file, height_map, x, y):
    """
    create 8 triangles from this pixel to the surrounding pixels
    
    triangles need to be wound counter clockwise
    """

    # compute the coordinates for reference later
    
    # scaled x - 1
    sxm1 = scale_pixel(x-1)
    # scaled x
    sx = scale_pixel(x)
    # scaled x + 1
    sxp1 = scale_pixel(x+1)
    
    # scaled y - 1
    sym1 = scale_pixel(y-1)
    # scaled y
    sy = scale_pixel(y)
    # scaled y + 1
    syp1 = scale_pixel(y+1)

    coords = {
        "n":    {"x": sx,    "y": syp1,   "z": height_map[x][y+1]},
        "ne":   {"x": sxp1,  "y": syp1,   "z": height_map[x+1][y+1]},
        "e":    {"x": sxp1,  "y": sy,     "z": height_map[x+1][y]},
        "se":   {"x": sxp1,  "y": sym1,   "z": height_map[x+1][y-1]},
        "s":    {"x": sx,    "y": sym1,   "z": height_map[x][y-1]},
        "sw":   {"x": sxm1,  "y": sym1,   "z": height_map[x-1][y-1]},
        "w":    {"x": sxm1,  "y": sy,     "z": height_map[x-1][y]},
        "nw":   {"x": sxm1,  "y": syp1,   "z": height_map[x-1][y+1]},
        "p":    {"x": sx,    "y": sy,     "z": height_map[x][y]},
    }


    # North pair
    # NNW
    facet(output_file, triangle(coords, "p", "n", "nw"))

    # NNE
    facet(output_file, triangle(coords, "p", "ne", "n"))
    
    # East pair
    # ENE
    facet(output_file, triangle(coords, "p", "e", "ne"))

    # ESE
    facet(output_file, triangle(coords, "p", "se", "e"))

    # South pair
    # SSW
    facet(output_file, triangle(coords, "p", "sw", "s"))

    # SSE
    facet(output_file, triangle(coords, "p", "s", "se"))

    # West pair
    # WSW
    facet(output_file, triangle(coords, "p", "w", "sw"))

    # WNW
    facet(output_file, triangle(coords, "p", "nw", "w"))
    
main()


