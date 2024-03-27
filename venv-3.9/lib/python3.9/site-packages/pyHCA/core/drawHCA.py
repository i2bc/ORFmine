#!/usr/bin/env python
""" The drawHCA module is an ensemble of functions to draw HCA plot of protein
sequences (with or without domain annotation
"""

import os, sys, argparse, string
import Bio
import Bio.SeqIO
import numpy as np

from pyHCA.core.seq_util import transform_seq, check_if_msa
from pyHCA.core.seq_util import compute_conserved_positions, compute_hca_conserved_positions
from pyHCA.core.ioHCA import read_annotation

__author__ = "Tristan Bitard-Feildel"
__licence__= "MIT"
__version__ = 0.1
__email__ = "tristan.bitard-feildel [you know what] impmc.upmc.fr"
__institute__ = "IMPMC"


def writeSVGheader(nbaa, fout):
    """
    brief write header of svg file
    param fout is an object file open in writing
    """
#<svg width="%d" height="4100" version="1.1" viewBox="0 0 %d 4"
    w = (((nbaa)/10.0)*90+150)*10
    fout.write("""<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg viewBox="0 0 60 55" width="200" height="100"
xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">

 <!-- ECMAScript with each click -->
  <script type="application/ecmascript"> <![CDATA[
    function aa_click(evt) {
      var aa = evt.target;
      var currentColor = aa.getAttribute("style");
      if (currentColor == "fill:blue; fill-opacity:0.0")
        {aa.setAttribute("style", "fill:blue; fill-opacity:0.3");
        }
      if (currentColor == "fill:blue; fill-opacity:0.3")
        {aa.setAttribute("style", "fill:red; fill-opacity:0.3");
        }
      if (currentColor == "fill:red; fill-opacity:0.3")
        {aa.setAttribute("style", "fill:gray; fill-opacity:0.3");
        }
      if (currentColor == "fill:gray; fill-opacity:0.3")
        {aa.setAttribute("style", "fill:blue; fill-opacity:0.0");
        }
        
    }
  ]]> </script>



""")
#%(w, int(w/100)))
    
def getSVGheader(nbaa, height=4100):
    """
    brief write header of svg file
    param fout is an object file open in writing
    """
    
    return ("""<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="{}" height="{}" version="1.1"
xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">

 <!-- ECMAScript with each click -->
  <script type="application/ecmascript"> <![CDATA[
    function aa_click(evt) {{
      var aa = evt.target;
      var currentColor = aa.getAttribute("style");
      if (currentColor == "fill:blue; fill-opacity:0.0")
        {{aa.setAttribute("style", "fill:blue; fill-opacity:0.3");
        }}
      if (currentColor == "fill:blue; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:red; fill-opacity:0.3");
        }}
      if (currentColor == "fill:red; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:gray; fill-opacity:0.3");
        }}
      if (currentColor == "fill:gray; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:blue; fill-opacity:0.0");
        }}
        
    }}
  ]]> </script>


""".format((((nbaa)/20.0)*90+150), height))

def getborder(x, y, F, yoffset=0):
    """
    """
    
    side  =    [[[ 0.83334, 7.555554,], [ -4.166666, 4.222221]],\
                [[-4.166666, 4.222221], [-2.5, -1.333335]],\
                [[-2.5, -1.333335], [4.166667, -3.555557]],\
                [[0.83334, 7.555554], [7.5, 5.333332]],\
                [[7.5, 5.333332],  [9.166667, -0.222224]],\
                [[9.166667, -0.222224], [4.166667, -3.555557]]] 
    
    lside = []
    polygon = "<polyline points=\""
    for begin, end in side:
        xb, yb, = begin
        xe, ye = end
        
        nxb = (x+xb)*F
        nyb = (y+yb*-1)*F
        nxe = (x+xe)*F
        nye = (y+ye*-1)*F
        
        polygon += "%f,%f %f,%f "%(nxb, nyb+yoffset, nxe, nye+yoffset)
        lside.append([[nxb,nyb], [nxe,nye]])
    polygon += """ " style="fill:black; fill-opacity:0.2" />"""
    return lside, polygon


def getCoordBorder(x, y, color="black", opacity=0.0, F=1, bis = 0, yoffset=0):
    """
    brief get the size limit of an amino acide in a HCAplot
    param x is the coord of the amino acid letter
    param y is the coord of the amino acid letter
    param F is the zoom index
    return lside is a list of 6 cases with coordinate of each x and y line around the amino acid in HCA plot
    """
    side  =    [[[ 0.83334, 7.555554,], [ -4.166666, 4.222221]],\
                [[-4.166666, 4.222221], [-2.5, -1.333335]],\
                [[-2.5, -1.333335], [4.166667, -3.555557]],\
                [[0.83334, 7.555554], [7.5, 5.333332]],\
                [[7.5, 5.333332],  [9.166667, -0.222224]],\
                [[9.166667, -0.222224], [4.166667, -3.555557]]] 
    
    
    lside = []
    polygon = "<polyline points=\""
    for begin, end in side:
        xb, yb, = begin
        xe, ye = end
        
        nxb = F * (x + xb)
        nyb = F * (y + yb) 
        nxe = F * (x + xe)
        nye = F * (y + ye) 
        lside.append([[nxb,nyb], [nxe,nye]])
        
        polygon += "%f,%f %f,%f "%(nxb+80, -nyb+80+bis+yoffset, nxe+80, -nye+80+bis+yoffset)
    
    polygon += """ " visibility="visible" style="fill:{}; fill-opacity:{}" onclick="aa_click(evt)"/>""".format(color, opacity)
    return lside, polygon
        
def getsideborder(seq, hydrophobe, nbaa):
    """
    brief get for each hydrophobe the cluster limit in a table
    param seq is the sequence with 4 spaces at the begining and at the end
    param hydrophobe is the string containing the hydrophobe amino acids
    param nbaa is the number of the initial protein
    return seqside is a table where each case contains 6 booleans representing the cluster limit
    """
    
    # adapted from the Luc Canard Algorithm
    
    # initialize of seqside to contains limit of each AA
    seqside = []
    for i in range(len(seq)):
        seqside.append([False, False, False, False, False, False])
    
    # for all position
    for n in range(0, len(seq)):
        if n >= 4 and n+1 < len(seq):       
            """ SIDE 0 """
            if seq[n] in hydrophobe:
                if seq[n-1] in hydrophobe:
                    seqside[n][0] = False
                else:
                    if seq[n-4] in hydrophobe:
                        if seq[n-3] == 'P' or seq[n-2] == 'P' or seq[n-1] == 'P':
                            seqside[n][0] = True
                        else:
                            seqside[n][0] = False
                    else:
                        seqside[n][0] = True
            else:
                if seq[n-4] in hydrophobe:
                    if seq[n-1] in hydrophobe:
                        if seq[n-3] == 'P' or seq[n-2] == 'P':
                            seqside[n][0] = False
                        else:
                            seqside[n][0] = True
                    else:
                        seqside[n][0] = False       
            """ SIDE 1 """
            if ( seq[ n] in hydrophobe):      
                if ( seq[ n - 4] in hydrophobe):
                    if seq[ n - 3] == 'P' or seq[ n - 2] == 'P' or seq[ n - 1] == 'P':
                        seqside[ n][1] = True
                    else:
                        seqside[ n][1]= False
                else:
                    if ( seq[ n - 3] in hydrophobe ):
                        if seq[ n - 2]  == 'P' or  seq[ n - 1] == 'P':
                            seqside[ n][1] = True
                        else:    
                            seqside[ n][1] = False
                    else:
                        seqside[ n][1] = True
            else:
                if (seq[n - 4] in hydrophobe):
                    if (seq[n - 3] in hydrophobe):
                        seqside[n][1] = True
                    else:
                        seqside[n][1] = False
                else:
                    seqside[n][1] = False
            
            """ SIDE 2"""
            if ( seq[ n] in hydrophobe):
                if ( seq[ n + 1] in hydrophobe):
                    seqside[ n][2] = False
                else:
                    if ( seq[ n - 3] in hydrophobe):
                        if ( seq[ n - 2] == 'P' or seq[ n - 1] == 'P'):    
                            seqside[ n][2] = True
                        else:
                            seqside[ n][2] = False
                    else:    
                        seqside[ n][2] = True
            else:
                if (seq[n - 3] in hydrophobe):
                    if (seq[n + 1] in hydrophobe):
                        if (seq[n - 2] == 'P' or seq[n - 1] == 'P'  or seq[n]    == 'P'    ):
                            seqside[n][2] = False
                        else:
                            seqside[n][2] = True
                    else:
                        seqside[n][2] = False
                else:
                    seqside[n][2] = False    

    
        if n < (len(seq)-4):
            """ SIDE 3"""
            if ( seq[ n] in hydrophobe ):
                if ( seq[ n - 1] in hydrophobe):
                    seqside[ n][3] = False
                else:
                    if (seq[ n + 3] in hydrophobe):
                        if ( seq[ n + 1]  == 'P' or seq[ n + 2] == 'P'):
                            seqside[ n][3] = True
                        else:
                            seqside[ n][3] = False
                    else:
                        seqside[ n][3] = True
            else:
                if (seq[n - 1] in hydrophobe):
                    if (seq[n + 3] in hydrophobe):
                        if (seq[n]  == 'P' or seq[n + 1]  == 'P' or  seq[n + 2] == 'P'):
                            seqside[n][3] = False
                        else:
                            seqside[n][3] = True
                    else:
                        seqside[n][3] = False
                
                else:
                    seqside[n][3] = False

            """ SIDE 4"""
            if ( seq[ n] in hydrophobe):
                if ( seq[ n +3] in hydrophobe):
                    if (     seq[ n + 1]  == 'P'   or seq[ n + 2] == 'P'):
                        seqside[n][4] = True
                    else:
                            seqside[n][4] = False
                else:
                    if ( seq[ n + 4] in hydrophobe):
                        if (    seq[ n + 1] == 'P' or seq[ n + 2] == 'P' or seq[ n + 3] == 'P'    ):
                            seqside[n][4] = True
                        else:
                            seqside[n][4] = False
                    else:
                        seqside[n][4] = True
            else:
                if n==1:
                    pass
                    #print n, seq[n], seq[n+3], seq[n+4]
                if ( seq[ n + 3] in hydrophobe):
                    if ( seq[ n + 4] in hydrophobe):
                        
                        seqside[n][4] = True
                        if n == 1 :
                            #print seqside[n]
                            pass
                    else:
                        seqside[n][4] = False
                else:
                    seqside[n][4] = False
        
            """ SIDE 5"""
            if ( seq[ n] in hydrophobe):
                if ( seq[ n + 1] in hydrophobe):
                    seqside[n][5] = False
                else:
                    if ( seq[ n + 4] in hydrophobe):
                        if ( seq[ n + 1] == 'P' or  seq[ n + 2] == 'P' or seq[ n + 3] == 'P' ):
                            seqside[n][5] = True
                        else:
                            seqside[n][5] = False
                    else:
                        seqside[n][5] = True
        
            else:
                if ( seq[ n + 4] in hydrophobe):
                    if ( seq[ n + 1] in hydrophobe):
                        if (    seq[ n + 2] == 'P'  or  seq[ n + 3] == 'P'):
                            seqside[n][5] = False
                        else:
                            seqside[n][5] = True
                    else:
                        seqside[n][5] = False
                else:
                    seqside[n][5] = False
            
    # I do not understand this part
    #around = [False] * 9
    #for n in range(nbaa+4, nbaa+2*4):
        
        #around[0] = (seq[n] in hydrophobe) and  (seq[n - 4] in hydrophobe)
        #around[1] = (seq[n]  in hydrophobe) and (seq[n - 3]  in hydrophobe)
        #around[3] = (seq[n]  in hydrophobe) and (seq[n - 1]  in hydrophobe)

        #for i in range(2, 6):
            #seqside[n][i] = False
        
        #if around[0]  and  around[3]:
            #seqside[n][0] = True
        #else:
            #seqside[n][0] = False
        #if around[0] and  around[3]:
            #seqside[n][1] = True
        #else:
            #seqside[n][1] = False

    #for n in range(0, 4):
        
        #around[5] = (seq[n] in hydrophobe) and (seq[n + 4] in hydrophobe)
        #around[7] = (seq[n] in hydrophobe) and (seq[n + 3] in hydrophobe)
        #around[8] = (seq[n] in hydrophobe) and (seq[n + 1] in hydrophobe)

        #for i in range(0, 4):
            #seqside[n][i] = False
        
        #if around[7]  and  around[8]:
            #seqside[n][ 4] = True
        #else:
            #seqside[n][ 4] = False
        #if around[5] and  around[8]:
            #seqside[n][ 5] = True
        #else:
            #seqside[n][ 5] = False

    return seqside

def drawside(side, dybis = 0, second=80, yoffset=0): #seqside, lside, n, seq):
    """
    brief draw line for cluster in HCAplot
    param side is list of coord x and y
    param dybis is for mirror
    """
    
    begin, end = side
    x1, y1 = begin
    x2, y2 = end
    
    x1 += 80
    x2 += 80
    
    y1 = -y1 +(dybis+second)
    y2 = -y2 + (dybis+second)
    
    return """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:blue;stroke:black;stroke-width:1;" />"""%(x1,y1+yoffset, x2,y2+yoffset)

def deplace(i, seq, coord, side, hydrophobe):
    """
    brief deplace or not the line cluster to avoid hole in line cluster
    param i is the iterator on seq
    param seq is the sequence with 4 spaces at the begining and at the end
    param coord is a list with n case containing x,y coordinate of amino acid in hca plot
    param side is the current side of line cluster
    param hydrophobe is a list containing list of aa hydrophobe
    return answer is the distance to deplace line cluster
    """
    
    """ Algorithm by Luc Canard transcription from C program """
    
    before = i - 1
    after = i + 1

    answer = 0.
    dy = 0.
    dymoins3 = 0.
    dyplus3 = 0.
    dyplus4 = 0.
    dymoins4 = 0.
    HEIGTH = 40
    if (seq[i] not in hydrophobe):
        dymoins4= coord[i][1] - coord[i - 4][1]
        dyplus3 = coord[i][1] - coord[i + 3][1] 
        dymoins3= coord[i][1] - coord[i - 3][1]
        dyplus4 = coord[i][1] - coord[i + 4][1]  
    
    #for j in range(6):
        if dymoins3 < 0.:
            dymoins3 = -dymoins3
        if dyplus3  < 0.:
            dyplus3 = -dyplus3
        if dymoins4 < 0.:
            dymoins4 = -dymoins4
        if dyplus4  < 0.:
            dyplus4 = -dyplus4
    
        if ((side == 0) and seq[before] in hydrophobe and (dymoins4 > 15.)):
            dy = coord[i][1] - coord[before][1]
            if (dy > 15.):
                answer = -HEIGTH
    
        if ((side == 3) and seq[before] in hydrophobe and (dyplus3 > 15.)):
            dy = coord[i][1] - coord[before][1]
            if (dy > 15.):
                answer = -HEIGTH
        
        if ((side == 2) and seq[after] in hydrophobe and (dyplus4 > 15.)):
            dy = coord[i][1] - coord[after][1]
            if (dy < -15.):
                answer = HEIGTH
        
        if ((side == 5) and seq[after] in hydrophobe and (dymoins3 > 15.)):
            dy = coord[i][1] - coord[after][1]
            if (dy < -15.):
                answer = HEIGTH
    
    return answer

def getCoord_pos(pos):
    """
    brief get coordinate of the plan helix for each position of the sequence 
    param seq is the sequence with 4 spaces at the end and at the begining
    return coord is a list with coordinate [x, y]
    """
    
    F = 1 # factor zoom
    FX = 12
    FY = 40
    NTOUR = 3.6
    BULGARIANCONSTANT = 11
    dx = FX / NTOUR
    dy = FY / NTOUR
        
    x = dx * (pos)                           # +80
    y = -dy * ((pos+BULGARIANCONSTANT)%NTOUR)  #    +80
        
    return x, y

def getCoord(seq, F=1, FX=12, FY=40):
    """
    brief get coordinate of the plan helix for each position of the sequence 
    param seq is the sequence with 4 spaces at the end and at the begining
    return coord is a list with coordinate [x, y]
    """
    
    AROUND = 4
    NTOUR = 3.6
    BULGARIANCONSTANT = 11
    dx = FX / NTOUR
    dy = FY / NTOUR
    #print(FX, FX/NTOUR,  12/NTOUR)
    coord = []
    #int_coord = []
    for n in range(len(seq)):
        
        x =  dx *  (n - AROUND)                             # +80
        y = -dy * ((n - AROUND + BULGARIANCONSTANT)%NTOUR)  # +80
        #ys = dy * ((n - AROUND + BULGARIANCONSTANT)%NTOUR) # +80 - 40
        #print(x, y, n-AROUND, ((n - AROUND + BULGARIANCONSTANT)%NTOUR))
        coord.append([x, y])
        #ix = int(round(x))
        #iy = int(round(y))
        #int_coord.append((ix, iy))
        
    coord = np.array(coord)
    return coord #, int_coord


### Sequence related drawing fonctions

def drawProline(x, y, yoffset=0):
    """
    brief draw a star for proline
    param x coordinate
    param y coordinate
    return polygone is a string to make svg star
    """
    
    displace = [[350,75],\
                [379,161],\
                [469,161],\
                [397,215],\
                [423,301],\
                [423,301],
                [350,250],\
                [277,301],\
                [303,215],
                [231,161],\
                [321,161]]
    
    polygon = "<polyline points=\""
    for cur_x, cur_y in displace:
        
        x1 = (cur_x-321)/22.0 + x
        y1 = (cur_y-250)/22.0 + y
        
        polygon += "%f,%f "%(x1, y1+yoffset)
    
    polygon += """ " style=" fill:black; fill-opacity:1" />"""
    
    return polygon  

def drawThreonine(x, y, yoffset=0):
    """
    brief draw a rectangle for threonine hca plot
    param x coord of x
    param y coord of y
    """
    
    rectangle = """<rect x="%f" y="%f" width="4" height="5" fill="none" stroke="black" stroke-width="0.2"/>"""%(x, y-5+yoffset)
    
    return rectangle
    
def drawSerine(x, y, yoffset=0):
    """
    brief draw a rectangle for serine hca plot
    param x coord of x
    param y coord of y
    """
    
    rectangle = """<rect x="%f" y="%f" width="4" height="5" fill="none" stroke="black" stroke-width="0.2"/>"""%(x, y-5+yoffset)
    prectangle = """<rect x="%f" y="%f" width="1" height="1" fill="black" stroke="black" stroke-width="0.2"/>"""%(x+1.5, y-3+yoffset)
    return rectangle+prectangle


def drawGlycine(x, y, yoffset=0):
    """
    brief draw a losange for glycine hca plot
    param x coord of x
    param y coord of y
    """
    displace = [[2,0],\
                [4,2],\
                [2,4],\
                [0,2]]
    
    polygon = "<polyline points=\""
    for cur_x, cur_y in displace:
        
        x1 = cur_x*1.2 + x+1
        y1 = cur_y*1.2 + y-5
        
        polygon += "%f,%f "%(x1, y1+yoffset)
    
    polygon += """ " style=" fill:black; fill-opacity:1" />"""
    return polygon
   

def linkCluster(seq, coords, n, dx, dy, yoffset=0):
    """
    brief trace the link to inside cluster
    param seq is a list containing residus and 4 spaces at the beg and at the end
    param coords is a list length of protein + 4 spaces in beg and 4 in end, each case contains [x, y] coordinates
    param n is the number of residu +4
    param dx to centre the draw
    param dy to centre the draw
    param fout svg opened file
    """
    
    dxmos = 4.0
    dymos = -3.5
    dxmos2 = 4.0
    dymos2 = 6.0
    
    svg = ""
    
    if coords[n+1][1] < coords[n][1]:
        x1 = coords[n][0]+dx+dxmos
        y1 = coords[n][1]+dy+dymos
        x2 = coords[n+1][0]+dx+dxmos
        y2 = coords[n+1][1]+dy+dymos
        svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)

        if coords[n+2][1] < coords[n+1][1]:
            x1 = coords[ n + 1][0] + dx + dxmos
            y1 = coords[ n + 1][1] + dy + dymos - 10
            x2 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos) * 0.5
            y2 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos - 40.) * 0.5
            svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        else:
            x1 = coords[ n + 1][0] + dx + dxmos
            y1 = coords[ n + 1][1] + dy + dymos + 40. - 10
            x2 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos) * 0.5
            y2 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos + 40.) * 0.5
            svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        
        x1 = coords[ n][0] + dx + dxmos
        y1 = coords[ n][1] + dy + dymos + 10
        x2 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos) * 0.5
        y2 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos + 40.) * 0.5
        svg +=  """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
    else:
        x1 = coords[ n][0] + dx + dxmos2
        y1 = coords[ n][1] + dy + dymos2
        x2 = coords[ n + 1][0] + dx + dxmos2
        y2 = coords[ n + 1][1] + dy + dymos2 - 20.
        svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        
        x1 = coords[ n][0] + dx + dxmos
        y1 = coords[ n][1] + dy + dymos
        x2 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos ) * 0.5
        y2 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos - 40.) * 0.5
        svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        
        x1 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos )  * 0.5
        y1 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos + 40.) * 0.5
        x2 = coords[ n + 1][0] + dx + dxmos
        y2 = coords[ n + 1][1] + dy + dymos
        svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        
    return svg
        

def dosvg(seq, coord, seqside, hydrophobe, conservation, nbaa, F=1, yoffset=0, idx_offset=0):
    """
    brief draw a svg hca plot
    param seq is a list containing sequence and 4 spaces at the begining and at the end
    param coord is a list  [x, y] amino acid coordinate for each position
    param seqside is a boolean list containing 6 cases by case. 6 faces True or False to draw line cluster
    param nbaa is the number of residue
    param idx_offset if the begining if we use cutting option
    """

    F = 1 # factor zoom
    FX = 12
    FY = 40
    AROUND = 4
    NTOUR = 3.6
    BULGARIANCONSTANT = 11
    dx = FX / NTOUR
    dy = FY / NTOUR

    coordbis = []
    for i, j in coord:
        coordbis.append([i, j-40])
    
    outsvg = ""
    for n in range(0, len(seq)):
        idx_seq = n - 4
        #print(idx_seq)
        color="black"
        opacity=0.0
        info_annot_position = conservation.get(idx_seq, dict())
        if "poly" in info_annot_position:
            color, opacity = info_annot_position["poly"]
        
        lside, polygon = getCoordBorder(coord[n][0], coord[n][1], color, opacity, yoffset=yoffset)
        lsides, polygon2 = getCoordBorder(coordbis[n][0], coordbis[n][1], yoffset=yoffset)
        
        x, y = coord[n]
        x += 80
        y = -y+80
        
        _, ys = coordbis[n]
        ys = -ys+80
        
        # draw letter
        if seq[n] == "P":
            outsvg += drawProline(x, y, yoffset)
            outsvg += drawProline(x, ys, yoffset)
        elif seq[n] == "S":
            outsvg += drawSerine(x, y, yoffset)
            outsvg += drawSerine(x, ys, yoffset)
        elif seq[n] == "T":
            outsvg += drawThreonine(x, y, yoffset)
            outsvg += drawThreonine(x, ys, yoffset)
        elif seq[n] == "G":
            outsvg += drawGlycine(x, y, yoffset)
            outsvg += drawGlycine(x, ys, yoffset)
        elif seq[n] == "C":
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman;stroke:black;stroke-width:0.4">%s</text>\n"""%(x*F, (y+yoffset)*F, seq[n])
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman;stroke:black;stroke-width:0.4">%s</text>\n"""%(x*F, (ys+yoffset)*F, seq[n])
        else:
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman">%s</text>\n"""%(x*F, (y+yoffset)*F, seq[n])
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman">%s</text>\n"""%(x*F, (ys+yoffset)*F, seq[n])
                
        # <ruler>
        if  n> 5 and ( n - 4 + 1) % ( 10) == 0 and n>=3:
            
            x =  coord[n][0] + 80 + 3.5#+80+dx
            dy1 = dy + coord[4][1] + 60 
            position = n - AROUND + 1 + idx_offset
            
            #print "AAAA",coord[n][0], dx
            #print "DIST", x-xp, x
            
            
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman">%s</text>\n"""%(x*F, dy1-5+yoffset, position)
            outsvg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x, dy1+yoffset, x, dy1+5+yoffset)
        
        # cluster linker
        if seq[n] in hydrophobe and n+3 < len(seq):
            if (seq[n+1] not in hydrophobe and seq[n+2] in hydrophobe and
               seq[n+3] not in hydrophobe and seq[n+1] != "P" and
               seq[n-1] not in hydrophobe):
               outsvg += linkCluster(seq, coord, n, 80, -80, yoffset)
               outsvg += linkCluster(seq, coordbis, n, 80, -80, yoffset)
                    
        # each side of amino acid        
        for iteside, is_side in enumerate(seqside[n]):
            if is_side == False:
                continue
            
            if n+4 >= len(seq):
                continue
            dybis = deplace(n, seq, coord, iteside, hydrophobe)
            dybis2 = deplace(n, seq, coordbis, iteside, hydrophobe)
            
            outsvg += drawside(lside[iteside], -dybis, yoffset=yoffset)
            outsvg += drawside(lsides[iteside], -dybis2, yoffset=yoffset)
        
        outsvg += polygon
        outsvg += polygon2
    return outsvg


### Domain related drawing fonctions

def domains2svg(start, stop, dom_name, status,  coord, yoffset=0):
    # ( n - 4 + 1) % ( 10) == 0 
    if start-1+4-1 > 0:
        x1 = coord[start+4-1][0] + 80 + 3.5#+80+dx
    else:
        x1 = coord[0][0] + 80 + 3.5#+80+dx
    if stop+4-1 < len(coord):
        x2 = coord[stop+4-1][0] + 80 + 3.5#+80+dx
    else :
        x2 = coord[-1][0] + 80 + 3.5
    size = x2 - x1
    
    name = """<text x="{}" y="{}" """.format(x1+10, yoffset+40)
    rect = """<rect x="{}" y="{}" """.format(x1, 25+yoffset)
    subrect = """<rect x="{}" y="{}" """.format(x1, 68+yoffset)
    if status == "!":
        name += """style="fill:black;font-size:12px;font-family:Times-Roman">"""
        rect += """style="fill:white; stroke:black; fill-opacity:1.0"  """
        subrect += """style="fill:white; stroke:black; fill-opacity:0.0"  """
    elif status == "o":
        name += """style="fill:black;font-size:12px;font-family:Times-Roman">"""
        rect += """style="fill:blue; stroke:black; fill-opacity:0.05"  """
        subrect += """style="fill:blue; stroke:black; fill-opacity:0.05"  """
    else:
        name += """style="fill:grey;font-size:12px;font-family:Times-Roman">"""
        rect += """style="fill:white; stroke:grey; stroke-dasharray=10,10; """
        subrect += """style="fill:white; stroke:grey; stroke-dasharray=10,10; """
        rect += """fill-opacity:1.0" """
        subrect += """fill-opacity:0.0" """
        
    name += """{}</text>\n""".format(dom_name)
    rect += """width="{}" height="20"/>\n\n""".format(size)
    subrect += """width="{}" height="100"/>\n\n""".format(size)
    return subrect+rect+name
    
#### computation of positions


def createHCAsvg(seq, nbaa, domains, conservation, yoffset=0, idx_offset=0):
    """
    brief draw the svg of hca
    param seq is the sequence with 4 spaces at both begining and end
    param nbaa is the true number of aa
    param output is the address of output file
    param idx_offset is the position to begin (for the ruler)
    """
    hydrophobe  = "YIMLFVW"
    seq = "    "+seq+"    "
    # n case with 6 sides in booleans
    seqside = getsideborder(seq, hydrophobe, nbaa)
    
    # get coord for amino acid in hca plot
    coord = getCoord(seq)
    
    
    # draw the svg
    svg = dosvg(seq, coord, seqside, hydrophobe, conservation, nbaa, yoffset=yoffset, idx_offset=idx_offset)
    if domains:
        svg += "\n".join([domains2svg(start, stop, name, status, coord, yoffset) 
                              for start, stop, name, nested, status in domains])
    return svg
    


def make_svg(prot, prev_seq, annot, yoffset=0, idx_offset=0):
    """ create svg for a given sequence
    """
    seq = transform_seq(prev_seq)
    nbaa = len(seq)
    svg = ""
    if prot != "":
        #yoffset = cnt*230
        svg += """<text x="0" y="{}" """.format(yoffset+25)
        svg += """style="fill:black;font-size:12px;font-family:Times-Roman">"""
        svg += """{}</text>\n""".format(prot)
    #else:
        #yoffset = cnt*120
    svg += createHCAsvg(seq, nbaa, annot.get("domains", list()), annot.get("positions", dict()), yoffset=yoffset, idx_offset=idx_offset)
    return svg, nbaa

def draw_columns_lines(columns_prot, columns_prev_prot, cnt):
    """ draw lines between selected columns of the two proteins
    """
    y_offset_prev = cnt*180
    y_offset_middle = y_offset_prev + (90)
    y_offset = y_offset_prev + 180
    F = 1
    svg = ""
    for col in columns_prot:
        idx_prot = columns_prot[col]
        idx_prev_prot = columns_prev_prot[col]
        if idx_prot > -1 and idx_prev_prot > -1:
            x_prev, y_prev = getCoord_pos(idx_prev_prot)
            x_prev += 80 + 0.5 # corresponds to left padding and letter centering 
            
            # middle position
            y_prev_middle = (y_offset_middle+110)*F
            
            # starting positions below first hca drawing
            x_prev, y_prev = x_prev*F, (y_offset_prev+170)*F 
            # position above current hca drawing
            x, y = getCoord_pos(idx_prot)
            x += 80 + 0.5 # corresponds to left padding and letter centering 
            x, y = x*F, (y_offset+50)*F
            if x < x_prev:
                # curved
                svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x_prev, y_prev, x, y_prev_middle)
                # straight
                svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x, y_prev_middle, x, y)
            else:
                # straight
                svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x_prev, y_prev, x_prev, y_prev_middle)
                # curved
                svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x_prev, y_prev_middle, x, y)
    return svg
    
def drawing(dfasta, annotation, pathout, window=-1, offset=list()):
    """ draw multiple fasta sequences
    """
    ext = os.path.splitext(pathout)[1]
    svg, svgheader = drawing_svg(dfasta, annotation, window, offset)
    if ext == ".svg":
        save_svg(pathout, svg, svgheader)
    elif ext == ".png":
        save_png(pathout, svg, svgheader)
    else:
        print("Warning, unknown save extension")
        pass
        drawing_plot(dfasta, annotation, pathout, window)

def drawing_svg(dfasta, annotation, window=-1, idx_offset=()):
    """ draw hca plot on a svg file
    """
    svg = ""
    max_aa = 0
    cnt = 0
    prev_prot = None
    yoffset = 0
    for i, (prot, prot_seq) in enumerate(dfasta.items()):
        if isinstance(prot_seq, Bio.SeqRecord.SeqRecord):
            prot_seq= str(prot_seq.seq)
        elif not isinstance(prot_seq, str):
            raise ValueError("Unknown amino acid sequence format {} for prot {} ".format(type(prot_seq), prot))
        if window != -1:
            # modulo sequence length
            cur_prot = prot
            for s in range(0, len(prot_seq), window):
                subseq = prot_seq[s:s+window+4] # +4 correspond to hca offset
                offset = s + idx_offset[i] - 1
                if s != 0:
                    prot = ""
                    yoffset += 120
                cur_svg, nbaa = make_svg(prot, subseq, annotation.get(cur_prot, dict()), yoffset, offset)
                svg += cur_svg
                if nbaa > max_aa:
                    max_aa = nbaa
                if s == 0:
                    if prev_prot != None and "columns" in annotation[cur_prot]:
                        print("warning cannot draw line conservation between protein if window is different of -1")
                    #svg += draw_columns_lines(annotation[prot]["columns"], annotation[prev_prot]["columns"], cnt-1)
                    prev_prot = prot
            yoffset += 180
        else:
            cur_svg, nbaa = make_svg(prot, prot_seq, annotation.get(prot, dict()), yoffset, idx_offset[i] - 1)
            svg += cur_svg
            if nbaa > max_aa:
                max_aa = nbaa
            if prev_prot != None and "columns" in annotation[prot]:
                svg += draw_columns_lines(annotation[prot]["columns"], annotation[prev_prot]["columns"], cnt-1)
            prev_prot = prot
            cnt += 1
            yoffset += 180
    
    # analys the new annotated domain, selective pressure from PAML
    #evolution_rate(pathnt, params.pathtree)
    svgheader = getSVGheader(max_aa, yoffset)
    return svg, svgheader

def save_svg(pathout, svg, svgheader):
    with open(pathout, "w") as fout:
        fout.write(svgheader)
        fout.write(svg)
        fout.write("</svg>")

def save_png(pathout, svg, svgheader):
    from cairosvg import svg2png

    svg_code = """{}{}\n</svg>
""".format(svgheader, svg)
    svg2png(bytestring=svg_code ,write_to=pathout)
            
            
def colorize_positions(msa, seq, conservation, method="identity"):
    """ colorize positions according to position
    """
    size = len(msa)
    positions = dict()
    #if method == "rainbow":
        #norm = matplotlib.colors.Normalize(vmin=0, vmax=size)
        #palette = plt.get_cmap("viridis")
        #pos = 0
        #for i in range(len(msa)):
            #if msa[i] != "-":
                ##print( i, pos, conservation[pos])
                #if conservation[pos] > 0.8:
                    #rgb = palette(norm(i))[:3]
                    #c = matplotlib.colors.rgb2hex(rgb)
                    #positions[pos] = {"poly": (c, 0.5)}
                #else:
                    #positions[pos] = {"poly": ("black", 0.0)}
                #pos += 1
    if "identity":
        #print(conservation)
        for i in range(len(seq)):
            if conservation[i] >= 0.9:
                positions[i] = {"poly": ("red", 0.5)}
            elif conservation[i] >= 0.66:
                positions[i] = {"poly": ("orange", 0.5)}
            elif conservation[i] >= 0.5:
                positions[i] = {"poly": ("yellow", 0.5)}
            else:
                positions[i] = {"poly": ("black", 0.0)}
    else:
        raise ValueError("Unknown colorization method")
    return positions

def select_columns(column_scores, msa2seq, threshold=0.8):
    segments = dict()
    for col in column_scores:
        if column_scores[col] >= threshold:
            idx, isgap = msa2seq[col]
            if not isgap:
                segments[col] = idx
            else:
                segments[col] = -1
    return segments

def get_params():
    """ get command line ArgumentParser
    """
    parser = argparse.ArgumentParser(prog="{} {}".format(os.path.basename(sys.argv[0]), "draw"))
    parser.add_argument("-i", action="store", dest="fastafile", help="the fasta file", required=True)
    parser.add_argument("-w", action="store", dest="window", type=int, help="sequence len before breaking the sequence to the next plot "
                        "(-1 the whole sequence are used, minimum size is 80)", default=-1)
    parser.add_argument("--offsets", action="store", dest="offsets", help="add seq index offsets to plot amino acid index", nargs="+", type=int, default=[-1])
    parser.add_argument("-d", action="store", dest="domain", help="[optionnal] provide domain annoation")
    parser.add_argument("-f", action="store", dest="domformat", help="the domain file format", choices=["pfam", "seghca"])
    #parser.add_argument("--color-msa", action="store", choices=["rainbow", "identity"], dest="msacolor", help="method to use to color a MSA", default="rainbow")
    parser.add_argument("--cons-msa", action="store", choices=["aa", "hca"], dest="msacons", help="method to use to compare sequences (aa, hca)")
    parser.add_argument("-o", action="store", dest="outputfile", help="output file in svg format", required=True)
    parser.add_argument("--verbose", action="store_true", dest="verbose", help="print information")
    params = parser.parse_args()
    
    if params.window > -1 and params.window < 80 and params.verbose:
        print("window parameter (-w) must either be superior to 80 or -1 to have the full sequence on one line")
    return params


def main():
    # params 
    params = get_params()
    
    from pyHCA.core.ioHCA import read_multifasta
    # read fasta file
    dfasta = read_multifasta(params.fastafile)
    if len(params.offsets) == 1 and params.offsets[0] == -1:
        offsets = [0] * len(dfasta)
        params.offsets = offsets[:]
        
    assert len(dfasta) == len(params.offsets), "different number of offsets {} found for {} sequences".format(len(dfasta), len(params.offsets))
    # are we using an msa ?
    dmsa = dict()
    is_an_msa = check_if_msa(dfasta)
    if is_an_msa:
        # if a msa is provided store sequence without gap character in a new dict
        # store msa sequence to get conserved positions 
        for rec in dfasta:
            seq = dfasta[rec]
            dmsa[rec] = seq
            if isinstance(seq, Bio.SeqRecord.SeqRecord):
                seq= str(seq.seq)
            dfasta[rec] = transform_seq(seq)
    if params.verbose:
        if is_an_msa :
            print("Reading {} as a MSA".format(params.fastafile))
        else:
            print("Reading {} as a list of fasta sequences".format(params.fastafile))
    if not is_an_msa and params.msacons is not None:
        print(params.msacons)
        print("Warning, the input sequence is not a MSA, the --cons-msa will be ignored", file=sys.stderr)
        
    # get conserved position
    if is_an_msa and params.msacons == "aa":
        msa_conserved_per_prot, msa_conserved_per_column, msa2seq = compute_conserved_positions(dfasta, dmsa)
    elif is_an_msa and params.msacons == "hca":
        msa_conserved_per_prot, msa_conserved_per_column, msa2seq = compute_hca_conserved_positions(dfasta, dmsa)
    
    # compute hca annotation  
    annotation = {}
    if params.domain:
        if params.verbose:
            print("Reading provided annotation")
        domains = read_annotation(params.domain, params.domformat)
        for prot in domains:
            annotation.setdefault(prot, dict())
            annotation[prot]["domains"] = domains[prot]
            
    # define msa color annotation
    if is_an_msa:
        if params.verbose:
            print("Getting colors for each positions")
        for prot in dfasta:
            annotation.setdefault(prot, dict())
            annotation[prot]["positions"] = colorize_positions(dmsa[prot], dfasta[prot], msa_conserved_per_prot[prot]) #, method=params.msacolor)
            annotation[prot]["columns"] = select_columns(msa_conserved_per_column, msa2seq.get(prot, dict()), threshold=1.0)
    else:
        # fil a dummy dictionary
        for prot in dfasta:
            annotation.setdefault(prot, dict())
            
    # draw
    ext = os.path.splitext(params.outputfile)[1]
    drawing(dfasta, annotation, params.outputfile, params.window, params.offsets)
        
    sys.exit(0)
    
    
if __name__ == "__main__":
    main()


