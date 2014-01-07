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
    convert_height_map(height_map(args.input_file))
       
 
       

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
    
    for x in range(0,x_dim):
        for y in range(y_dim -1,0,-1):
            (r,g,b) = pix[x,y]
            
            # r+g+b = meters * scale factor(???)
            values[x][y_dim-y] = (r*0.1)+(g*0.01)+(b*-0.1)

    return values


def stl_header():
    print("solid DEM")
    
def stl_footer():
    print("endsolid")


def convert_height_map(height_map):
    """
    Convert height map to STL data


    STL axis order:


    ^ 
    |Y+
    |
    0,0,0 -------> x+


    ( y axis flipped vs jpeg)
    """

    stl_header()    

    for x in range(0,len(height_map)-1):
        for y in range(0,len(height_map[x])-1):

            # draw a pair of triangles (square) for each pixel equivalent
            print("facet normal 1 1 1")
            print("\touter loop")
            print("\t\tvertex %d %d %d" % (x,y, height_map[x][y]))
            print("\t\tvertex %d %d %d" % (x+1,y+1, height_map[x+1][y+1]))
            print("\t\tvertex %d %d %d" % (x,y+1, height_map[x][y+1]))
            print("\tendloop")
            print("endfacet")
                        
            print("facet normal 1 1 1")
            print("\touter loop")
            print("\t\tvertex %d %d %d" % (x,y, height_map[x][y]))
            print("\t\tvertex %d %d %d" % (x+1,y, height_map[x+1][y]))
            print("\t\tvertex %d %d %d" % (x+1,y+1, height_map[x+1][y+1]))
            print("\tendloop")
            print("endfacet")


    stl_footer()

main()


