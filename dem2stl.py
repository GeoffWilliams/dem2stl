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



def main():
    logging.info("DEM TO STL")
    parser = argparse.ArgumentParser(description='DEM (Digital Elevation Model) to STL (3D printer) conversion')
    parser.add_argument('--input-file', dest='input_file', action='store',
            help='filename for input file', required=True)
    parser.add_argument('--output-file', dest='output_file', action='store',
            help='filename for output file', required=True)
    parser.add_argument('--r-ratio', dest='r_ratio', action='store',
            help='ratio for red pixels', default=0.05, type=float)
    parser.add_argument('--g-ratio', dest='g_ratio', action='store',
            help='ratio for green pixels', default=0.01, type=float)
    parser.add_argument('--b-ratio', dest='b_ratio', action='store',
            help='ratio for blue pixels', default=-0.05, type=float)
    parser.add_argument('--pixel-ratio', dest='pixel_ratio', action='store',
            help='pixel to world unit ratio', default=1, type=float)
    parser.add_argument('--resample-size', dest='resample_size', action='store',
            help='resample the source image by averaging squares of this size (px)', type=int, default=1)
    parser.add_argument('--debug', dest='debug', action='store_true',
            help='enable debug messages', default=False)

    args = parser.parse_args()
    
    if (args.debug):
        print("enabling debug mode...")
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)
        logging.debug("...debug mode enabled")

    logging.debug("RGB and pixel scale ratios:  r=%f g=%f b=%f p=%f" % 
        (args.r_ratio, args.g_ratio, args.b_ratio, args.pixel_ratio)) 

    convert_height_map(
        height_map(args.input_file, args.r_ratio, args.g_ratio, args.b_ratio, args.resample_size), 
        args.output_file, args.pixel_ratio)
       
 
       

def height_map(input_file, r_ratio, g_ratio, b_ratio, resample_size):
    """
    Convert the image pixel values to a 2D array of heights



    pixel order JPEG:

    0,0     --------->X+
    |
    |
    | y+
    v
    """
    im = Image.open(input_file) #Can be many different formats.
    pix = im.load()

    x_dim, y_dim = im.size
    logging.debug("creating map %dx%d" % (x_dim, y_dim))
    values = numpy.empty((x_dim, y_dim))
 
    for x in range(0,x_dim):
        for y in range(y_dim -1,0,-1):
            (r,g,b) = pix[x,y]
            
            # r+g+b = meters * scale factor(???)
            values[x][y_dim-y] = (r*r_ratio)+(g*g_ratio)+(b*-b_ratio)

    if (resample_size > 1):
        values = resample(values, resample_size)

    return values

def average(values, resample_size, source_x_dim, source_y_dim, resampled_x, resampled_y):
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
    
    resample_squared = math.pow(resample_size * 2 + 1,2)
    logging.debug("resample size %d will average squares of %d pixels" % (resample_size, resample_squared))

    source_x = resampled_x * resample_size
    source_y = resampled_y * resample_size

    # find the x,y values for each corner of the square illustrated above
    x_min = source_x - resample_size
    x_max = source_x + resample_size
    y_min = source_y - resample_size
    y_max = source_y + resample_size
    
    # take the mean of the pixels indicated by resample size
    v=0
    for x in range(x_min, x_max):
        for y in range(y_min, y_max):
            v += values[x][y]     

    v = v / resample_squared
    return v


def resample(values, resample_size):
    source_x_dim = len(values)
    source_y_dim = len(values[0])

    # calculate resampled image size.  Round down (throw away) any pixels that
    # dont fit neatly
    x_dim = int(round(source_x_dim / resample_size))
    y_dim = int(round(source_y_dim / resample_size))

    logging.debug("resample to %f x %f" % (x_dim, y_dim))

    px_count = resample_size * resample_size
    resampled = numpy.empty((x_dim, y_dim))
    for x in range(0, x_dim):
        for y in range(0, y_dim):
            resampled[x][y]=average(values, resample_size, source_x_dim, source_y_dim, x, y)

    return resampled

def stl_header(output_file):
    output_file.write("solid DEM\n")
    
def stl_footer(output_file):
    output_file.write("endsolid\n")


def convert_height_map(height_map, output_filename, pixel_ratio):
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
            link_pixel(output_file, height_map, pixel_ratio, x, y)

    # link the N, S, E, W sides to zero height (the previously unprocessed border)
    
    # N
    for x in range(1, len(height_map)):
        link_floor_x(output_file, height_map, x)
    

    # S

    # E

    # W

    # link the bottom sides to create a solid
    

    stl_footer(output_file)

def link_floor_x(output_file, height_map, x):

    coords = {
        "n":    {"x": x,    "y": y+1,   "z": height_map[x][y+1]},
        "ne":   {"x": x+1,  "y": y+1,   "z": height_map[x+1][y+1]},
        "e":    {"x": x+1,  "y": y,     "z": height_map[x+1][y]},
        "se":   {"x": x+1,  "y": y-1,   "z": height_map[x+1][y-1]},
        "s":    {"x": x,    "y": y-1,   "z": height_map[x][y-1]},
        "sw":   {"x": x-1,  "y": y-1,   "z": height_map[x-1][y-1]},
        "w":    {"x": x-1,  "y": y,     "z": height_map[x-1][y]},
        "nw":   {"x": x-1,  "y": y+1,   "z": height_map[x-1][y+1]},
        "p":    {"x": x,    "y": y,     "z": height_map[x][y]},
    }


    # single triangle

    # double triangles
#    for (i in range(1,len(coords) - 1):
        
    # single triangle
    

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

def scale_pixel(pixel_ratio, n):
    return pixel_ratio * n


def link_pixel(output_file, height_map, pixel_ratio, x, y):
    """
    create 8 triangles from this pixel to the surrounding pixels
    
    triangles need to be wound counter clockwise
    """

    # compute the coordinates for reference later
    
    # scaled x - 1
    sxm1 = scale_pixel(pixel_ratio, x-1)
    # scaled x
    sx = scale_pixel(pixel_ratio, x)
    # scaled x + 1
    sxp1 = scale_pixel(pixel_ratio, x+1)
    
    # scaled y - 1
    sym1 = scale_pixel(pixel_ratio, y-1)
    # scaled y
    sy = scale_pixel(pixel_ratio, y)
    # scaled y + 1
    syp1 = scale_pixel(pixel_ratio, y+1)

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


