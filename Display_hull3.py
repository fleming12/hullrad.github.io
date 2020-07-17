"""
Display_hull3 (Uses Python Version 3.x for newer Python3 based Pymol versions)

PyMOL script to dispaly a solid convex hull of a selection of atoms.
AUTHOR: Patrick Fleming, 2019

Adapted from PyDeT for visualizing a Delaunay triangulation of a protein.
Reference:
R. Ordog PyDeT, a PyMOL plug-in for visualizing geometric concepts around proteins. Bioinformation
2(8): 346-347 (2008)

If your python3.x installation includes numpy and scipy YOU DON'T NEED ANYTHING ELSE.
This is the preferred way to go.

If you don't have numpy and scipy installed, you will need a separate program to
    calculate the convex hull - called qconvex.
This is a program in the qhull suite of programs.

Qhull may be downloaded from http://www.qhull.org/download.
        For UNIX download "Qhull_2015.2 for Unix"
        tar xzvf qhull-2015-src-7.2.0.tgz
        cd qhull-2015.2/
        make
        sudo make install
        and the path for qconvex will be /usr/local/bin/qconvex

OS X qhull binaries are available from MacPorts package manager.
        type "sudo port install qhull".
        and the path for qconvex will be /opt/local/bin/qconvex

OS X binaries are also available from Fink package manager.

The default in this script is /opt/local/bin/qconvex but you may
have to change it if you are using UNIX

After installing qhull go to: ### Edit the path to qconvex ### below (~line 120) and
change /opt/local/bin/qconvex to /usr/local/bin/qconvex if that is where qconvex is.

If you are using OS X but are familiar with UNIX and have gcc installed, the UNIX
instructions above work fine on a Mac.

USAGE: Deposit this script in the working directory, then
       in the PyMOL command line type, 
        run Display_hull.py 
        and follow the prompt:

        "hull [selection]" for hull with planes only
          or
        "hull [selection], True" for hull and edges as lines

Change the "ALPHA" variable below to change the transparency of the convex hull.

4/14/2020 SBP - Edits to get to work in PyMol 2.3.4 on MacOS.
(Many thanks to Shae Padrick for initial porting to Python Version 3)
"""

import sys
from pymol import cmd
from pymol.cgo import *
import string
import math
import shutil

#If useNumpy stays False, then it will look for qconvex
useNumpy = False
try:
    import numpy as np
    from scipy.spatial import ConvexHull
    useNumpy = True
except:
    pass

# Prompt to user
print('Enter "hull [selection],[True] to display a convex hull around selection +/- solid edges"')

def compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    # A helper function for computing the normal to a triangular facet
    nx = (y2-y1)*(z3-z2) - (z2-z1)*(y3-y2)
    ny = (z2-z1)*(x3-x2) - (x2-x1)*(z3-z2)
    nz = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)

    return (nx,ny,nz)

def draw_plane_cgo(obj, name,apex1,apex2,apex3,apex4,color=(0.5,0.5,0.5)):
  """
    Create a CGO plane from three arbitary coordinates

    Usage:
    draw_plane_cgo apex1, apex2, apex3, apex4, color

    where each apex is a 3-element vector and color is a 3-element RGB
    list defining the color of the plane (where each value of R, G
    and B is between 0 and 1 inclusive).
  """

  # Convert args to floating point numbers
  x1,y1,z1 = map(float,apex1)
  x2,y2,z2 = map(float,apex2)
  x3,y3,z3 = map(float,apex3)
  x4,y4,z4 = map(float,apex4)
  if type(color) == type(''):
    color = map(float,color.replace('(','').replace(')','').split(','))

  # Compute the normal vector for the triangle
  normal1 = compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3)
  normal2 = compute_normal(x1, y1, z1, x3, y3, z3, x4, y4, z4)
  normal3 = compute_normal(x2, y2, z2, x3, y3, z3, x4, y4, z4)

  # Create the CGO objects
  # Uncomment ALPHA line below to draw semi-transparent hull
  obj.extend([
    BEGIN, TRIANGLE_STRIP,
    COLOR, color[0], color[1], color[2],
    ALPHA, 0.3,
    NORMAL, normal1[0], normal1[1], normal1[2],
    VERTEX, x1, y1, z1,
    VERTEX, x2, y2, z2,
    VERTEX, x3, y3, z3,
    VERTEX, x4, y4, z4,

    END
  ])

def Hull(selection='all',edges = False,name=''):
    infile="pydelaunay.input"
    outfile="pydelaunay.output"
    statfile="pydelaunay.stat"

    ### Edit the path to qconvex ###
    qconvex= '/opt/local/bin/qconvex' #For Linux enviorment

    if useNumpy:
        Convex_SolidVisualise(outfile,selection,name)
        if edges:
            Convex_TessVisualise(outfile,selection,name)
    else:
        Convex_DExtract(infile,selection)
        Convex_QConvex(qconvex,infile,outfile,statfile)
        Convex_SolidVisualise(outfile,selection,name)
        if edges:
            Convex_TessVisualise(outfile,selection,name)

        print(' Returning statistics ')

        fpin = open(statfile)
        print('')
        print(fpin.read())

def Convex_DExtract(infile, selection='all'):
    print(' Selection in use: '+selection)
    print(' Extracting data for qconvex into file: '+infile)
    print(' -> Collecting data')
    model = cmd.get_model(selection)
    li='3\n'
    li=li+'%s\n'%cmd.count_atoms(selection)
    
    for a in model.atom:
        li=li+'%s %s %s\n'%(a.coord[0],a.coord[1],a.coord[2])
    
    print(' -> Saveing data')
    fpout = open(infile, 'w')
    fpout.write(li)
    fpout.close()
    
    print(' -> Data extracted')
    
def Convex_QConvex(qconvex,infile,outfile,statfile='', options='p'):
    if statfile=='':
        print(' Starting qconvex! Geometry output is: '+outfile)
        os.system(qconvex+' -i '+options+' < '+infile+' > '+outfile)
    else:
        print(' Starting qconvex! Geometry output is: '+outfile+' and stat output is: '+statfile)
        os.system(qconvex+' -s -i '+options+' < '+infile+' > '+outfile+' 2> '+statfile)
    
def Convex_SolidVisualise(outfile,selection,name=''):
    print(' Visualise the Convex Hull')
    
    model = cmd.get_model(selection)
    coords = []
    for atm in model.atom:
        coords.append((float(atm.coord[0]), float(atm.coord[1]), float(atm.coord[2])))
    if useNumpy:
        coords = np.array(coords)
    else:
        pass
    if useNumpy:
        obj = []
        name = ''
        model_np = np.array(model)
        convex_hull = ConvexHull(coords)
        for r1 in convex_hull.simplices:
            a_atom = r1[0]
            b_atom = r1[1]
            c_atom = r1[2]
            coor1 = coords[a_atom]
            coor2 = coords[b_atom]
            coor3 = coords[c_atom]
            coor4 = coords[a_atom]
            draw_plane_cgo(obj, name,coor1,coor2,coor3,coor4)

    else:
        fpin = open(outfile)
        lines = open(outfile).readlines()
        numlines = int(lines[0])
        print('    Number of regions: ',numlines)
        
        obj = []
        name = ''
        for i in range(1,numlines):
            line=lines[i].strip()
            r1=line.split(' ')
            a_atom = eval(r1[0])
            b_atom = eval(r1[1])
            c_atom = eval(r1[2])
            #coor1 = model.atom[a_atom].coord
            #coor2 = model.atom[b_atom].coord
            #coor3 = model.atom[c_atom].coord
            #coor4 = model.atom[a_atom].coord
            coor1 = coords[a_atom]
            coor2 = coords[b_atom]
            coor3 = coords[c_atom]
            coor4 = coords[a_atom]
            draw_plane_cgo(obj, name,coor1,coor2,coor3,coor4)
        
        fpin.close()

    cmd.bg_color('white')
    cmd.load_cgo(obj,'Hull')
    
def Convex_TessVisualise(outfile,selection,name=''):
    print(' Visualise the Convex Tessalation (Edges)')

    if useNumpy:
        model = cmd.get_model(selection)
        coords = []
        for atm in model.atom:
            coords.append((float(atm.coord[0]), float(atm.coord[1]), float(atm.coord[2])))
        if useNumpy:
            coords = np.array(coords)
        convex_hull = ConvexHull(coords)
        numlines = len(convex_hull.simplices)

        print('    Number of regions: ', numlines)
           
        edgeblocker=len(model.atom)*[0]
        for i in range(len(edgeblocker)):
            edgeblocker[i]=len(model.atom)*[0]
    
        edgelist=[]
            
        a = model.atom[0]
        b = model.atom[1]
        minco=math.sqrt((a.coord[0]-b.coord[0])**2+(a.coord[1]-b.coord[1])**2+(a.coord[2]-b.coord[2])**2)
        print('minco = ', minco)
        maxco=0
        for r1 in convex_hull.simplices:
            for i in range(len(r1)):
                for j in range(len(r1)):
                    ii=(r1[i])
                    jj=(r1[j])
                    if ii<jj:        
                        if edgeblocker[ii][jj]==0:
                            edgeblocker[ii][jj]=1
                            edgeblocker[jj][ii]=1
            
                            a = model.atom[ii]
                            b = model.atom[jj]
                            co=math.sqrt((a.coord[0]-b.coord[0])**2+(a.coord[1]-b.coord[1])**2+(a.coord[2]-b.coord[2])**2)
                            minco=min(co,minco)
                            maxco=max(co,maxco)                     
                            edgelist.extend([[a.coord[0],a.coord[1],a.coord[2],b.coord[0],b.coord[1],b.coord[2],co]])
            
    else:
        print(' -> Reading file '+outfile+' and Calculating distances')
    
        model = cmd.get_model(selection)
        fpin = open(outfile)
        lines = open(outfile).readlines()
        numlines = int(lines[0])
        print('    Number of regions: ', numlines)
        
        edgeblocker=len(model.atom)*[0]
        for i in range(len(edgeblocker)):
            edgeblocker[i]=len(model.atom)*[0]
    
        edgelist=[]
        
        a = model.atom[0]
        b = model.atom[1]
        minco=math.sqrt((a.coord[0]-b.coord[0])**2+(a.coord[1]-b.coord[1])**2+(a.coord[2]-b.coord[2])**2)
        print('minco = ', minco)
        maxco=0
        for i in range(1,numlines):
            line=lines[i].strip()
            corner=line.split(' ')
    
            for i in range(len(corner)):
                for j in range(len(corner)):
                    ii=eval(corner[i])
                    jj=eval(corner[j])
                    if ii<jj:
                        if edgeblocker[ii][jj]==0:
                            edgeblocker[ii][jj]=1
                            edgeblocker[jj][ii]=1
    
                            a = model.atom[ii]
                            b = model.atom[jj]
                            co=math.sqrt((a.coord[0]-b.coord[0])**2+(a.coord[1]-b.coord[1])**2+(a.coord[2]-b.coord[2])**2)
                            minco=min(co,minco)
                            maxco=max(co,maxco)
                            edgelist.extend([[a.coord[0],a.coord[1],a.coord[2],b.coord[0],b.coord[1],b.coord[2],co]])
    
        fpin.close()
    

    print(' -> Generating CGO list')
    difco=maxco-minco
    obj=[]
    for e in edgelist:
        #co=((e[6]-minco)/difco)**(0.75)
        
        #co_r=max(min(1-2*co,1),0)
        #co_g=max(min(1-abs(2*co-1),1),0)
        #co_b=max(min(2*co-1,1),0)
        co_r = 0.00
        co_g = 0.00
        co_b = 0.00
        obj.extend([ CYLINDER, e[0],e[1], e[2], e[3],e[4], e[5],0.05,co_r,co_g,co_b,co_r,co_g,co_b ])
    print(' -> Create CGO objects ')
    cmd.load_cgo(obj,'Edges')

cmd.extend("hull",Hull)

