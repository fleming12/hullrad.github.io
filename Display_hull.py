"""
Display_hull
PyMOL script to dispaly a solid convex hull of a protein structural model.

Adapted from PyDeT for visualizing a Delaunay triangulation of a protein.
Reference:
R. Ordog PyDeT, a PyMOL plug-in for visualizing geometric concepts around proteins. Bioinformation
2(8): 346-347 (2008)

Requires qconvex in path. This is a program in the qhull suite of programs.
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

The default in this script is the UNIX path (/usr/local/bin/qconvex) but you may
have to change it if you are using OS X.

After installing qhull go to: ### Edit the path to qconvex ### below (~line 107) and
change /usr/local/bin/qconvex to /opt/local/bin/qconvex if that is where qconvex is.

If you are using OS X but are familiar with UNIX and have gcc installed, the UNIX
instructions above work fine on a Mac.

USAGE: Deposit this script in the working directory, then
       in the PyMOL command line type, 
        run Display_hull.py 
        and follow the prompt.

Note that the convex hull will be a compiled graphics object in PyMOL and you will not be 
able to ray trace it.

Change the "ALPHA" variable below to change the transparency of the convex hull.
"""

from pymol import cmd
from pymol.cgo import *
import string
import math
import shutil

# Prompt to user
print 'Enter "shull [selection] to display a solid convex hull around selection"'

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
#   ALPHA, 0.5,
    NORMAL, normal1[0], normal1[1], normal1[2],
    VERTEX, x1, y1, z1,
    VERTEX, x2, y2, z2,
    VERTEX, x3, y3, z3,
    VERTEX, x4, y4, z4,

    END
  ])

def Hull(selection='all',name=''):
	infile="pydelaunay.input"
	outfile="pydelaunay.output"
	statfile="pydelaunay.stat"

	### Edit the path to qconvex ###
	qconvex= '/usr/local/bin/qconvex' #For Linux enviorment

	Convex_DExtract(infile,selection)
	Convex_QConvex(qconvex,infile,outfile,statfile)
	Convex_SolidVisualise(outfile,selection,name)

	print ' Returning statistics '

	fpin = open(statfile)
	print '' 
	print fpin.read()

def Convex_DExtract(infile, selection='all'):
	print ' Selection in use: '+selection
	print ' Extracting data for qconvex into file: '+infile
	print ' -> Collecting data'
	model = cmd.get_model(selection)
	li='3\n'
	li=li+'%s\n'%cmd.count_atoms(selection)
	
	for a in model.atom:
	   li=li+'%s %s %s\n'%(a.coord[0],a.coord[1],a.coord[2])
	
	print ' -> Saveing data'
	fpout = open(infile, 'w')
	fpout.write(li)
	fpout.close()
	
	print ' -> Data extracted'
	
def Convex_QConvex(qconvex,infile,outfile,statfile='', options='p'):
	if statfile=='':
	   print ' Starting qconvex! Geometry output is: '+outfile
	   os.system(qconvex+' -i'+options+' < '+infile+' > '+outfile)
	else:
	   print ' Starting qconvex! Geometry output is: '+outfile+' and stat output is: '+statfile
	   os.system(qconvex+' -s -i'+options+' < '+infile+' > '+outfile+' 2> '+statfile)
	
def Convex_SolidVisualise(outfile,selection,name=''):
	print ' Visualise the Convex Tessalation (Edges)'
	print ' -> Reading file '+outfile+' and Calculating distances'
	
	model = cmd.get_model(selection)
	fpin = open(outfile)
	print '    Number of regions: '+fpin.readline()+''
	
	obj = []
        name = ''
	for line in fpin.readlines():
		line=line.strip()
		corner=line.split(' ')
		a_atom = eval(corner[0])
		b_atom = eval(corner[1])
		c_atom = eval(corner[2])
        	coor1 = model.atom[a_atom].coord
        	coor2 = model.atom[b_atom].coord
        	coor3 = model.atom[c_atom].coord
        	coor4 = model.atom[a_atom].coord
		draw_plane_cgo(obj, name,coor1,coor2,coor3,coor4)
	
	fpin.close()
        basedir = os.getcwd()

	cmd.load_cgo(obj,name)
	
cmd.extend("shull",Hull)

