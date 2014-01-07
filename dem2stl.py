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
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import argparse
import Image
import numpy

def main():
    print("DEM TO STL")
    parser = argparse.ArgumentParser(description='DEM (Digital Elevation Model) to STL (3D printer) conversion')
    parser.add_argument('--input-file', dest='input_file', action='store',
            help='filename for input file', required=True)
    parser.add_argument('--output-file', dest='output_file', action='store',
            help='filename for output file', required=True)

    args = parser.parse_args()
    convert_height_map(height_map(args.input_file), args.output_file)
       
 
       

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
    im = Image.open(input_file) #Can be many different formats.
    pix = im.load()
#    print im.size #Get the width and hight of the image for iterating over

    x_dim, y_dim = im.size
#    print ("creating map %dx%d" % (x_dim, y_dim))
    values = numpy.empty((x_dim, y_dim))
 
    x_dim=3
    y_dim=3
   
    for x in range(0,x_dim):
        for y in range(y_dim -1,0,-1):
            (r,g,b) = pix[x,y]
            
            # r+g+b = meters * scale factor(???)
            values[x][y_dim-y] = (r*0.1)+(g*0.01)+(b*-0.1)

    return values


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

    for x in range(1,3):
#len(height_map) -1):
        for y in range(1,3):
#len(height_map[x]) -1):
            link_pixel(output_file, height_map, x, y)

    stl_footer(output_file)

def facet(output_file, triangle, normal="1 1 1"):
    """
    write a facet using arrays of coordinates
    """
    output_file.write("facet normal %s\n" % (normal))
    output_file.write("\touter loop\n")
    output_file.write("\t\tvertex %d %d %d\n" % (triangle["xs"][0], triangle["ys"][0], triangle["zs"][0]))
    output_file.write("\t\tvertex %d %d %d\n" % (triangle["xs"][1], triangle["ys"][1], triangle["zs"][1]))
    output_file.write("\t\tvertex %d %d %d\n" % (triangle["xs"][2], triangle["ys"][2], triangle["zs"][2]))
    output_file.write("\tendloop\n")
    output_file.write("endfacet\n")


def triangle(coords, k1, k2, k3):
    """
    construct a map of ordered coordinates
    """
    return {
        "xs": [coords[k1]["x"], coords[k2]["x"], coords[k3]["x"]], 
        "ys": [coords[k2]["y"], coords[k2]["y"], coords[k3]["y"]],
        "zs": [coords[k3]["z"], coords[k2]["z"], coords[k3]["z"]],
    }




def link_pixel(output_file, height_map, x, y):
    """
    create 8 triangles from this pixel to the surrounding pixels
    """

    # compute the coordinates for reference later
    coords = {
        "n":    {"x": x,    "y": y+1,   "z": height_map[x][y+1]},
        "ne":   {"x": x+1,  "y": y+1,   "z": height_map[x+1][y+1]},
        "e":    {"x": x+1,  "y": y,     "z": height_map[x+1][y]},
        "se":   {"x": x-1,  "y": y+1,   "z": height_map[x-1][y+1]},
        "s":    {"x": x,    "y": y-1,   "z": height_map[x][y-1]},
        "sw":   {"x": x-1,  "y": y-1,   "z": height_map[x-1][y-1]},
        "w":    {"x": x-1,  "y": y,     "z": height_map[x-1][y]},
        "nw":   {"x": x-1,  "y": y+1,   "z": height_map[x-1][y+1]},
        "p":    {"x": x,    "y": y,     "z": height_map[x][y]},
    }


    # North pair
    # NNW
    facet(output_file, triangle(coords, "p", "nw", "n"))

    # NNE
    facet(output_file, triangle(coords, "p", "n", "ne"))
    
    # East pair
    # ENE
    facet(output_file, triangle(coords, "p", "ne", "e"))

    # ESE
    facet(output_file, triangle(coords, "p", "e", "se"))

    # South pair
    # SSW
    facet(output_file, triangle(coords, "p", "se", "s"))

    # SSE
    facet(output_file, triangle(coords, "p", "s", "sw"))

    # West pair
    # WSW
    facet(output_file, triangle(coords, "p", "sw", "w"))

    # WNW
    facet(output_file, triangle(coords, "p", "w", "nw"))
    
main()


