import sys, os, subprocess
from sys import argv
from math import pi, cos, sin, sqrt, floor, ceil

# import all user-defined functions for working with water files
#import clusters
#import clusterdata
#import recursive

#========= constants ============== 
ang2bohr = 1.889725989
hartree2kjmol=2625.49962
kjmol2hartree=1.0/hartree2kjmol
# == TIP3P model of water
rOH = 0.9572 #ang
aHOH = 104.52 * pi / 180.0
charge = [-0.834, +0.417, +0.417]
sigma = 3.15061 # ang
epsilon = 0.6364 #kJ/mol

def GenerateCoords(params):
 
 coords = []

 R = params[0]
 Th = params[1]
 Fi = params[2]

 rOHxy = rOH * cos( aHOH/2.0 )
 rOHz  = rOH * sin( aHOH/2.0 )
 
 for iatom in range(0,6):

	xyz = [0.,0.,0.]
	if ( iatom == 0 ):
	 xyz = [0.,0.,0.]

	elif ( iatom == 1 ):
	 x = rOHxy * cos( Fi )
	 y = rOHxy * sin( Fi )
	 z = rOHz 
	 xyz = [x,y,z]

	elif ( iatom == 2 ):
	 x = rOHxy * cos( Fi )
	 y = rOHxy * sin( Fi )
	 z = -rOHz 
	 xyz = [x,y,z]

	elif ( iatom == 3 ):
	 xyz = [R,0.,0.]

	elif ( iatom == 4 ):
	 ang = Th - aHOH/2.0
	 x = R + rOH * cos( ang )
	 y = rOH * sin( ang )
	 xyz = [x,y,0.]

	elif ( iatom == 5 ):
	 ang = Th + aHOH/2.0
	 x = R + rOH * cos( ang )
	 y = rOH * sin( ang )
	 xyz = [x,y,0.]

	else:
	 print "Atoms number is too large"
	 exit(3)
	 
	coords.append(xyz[:])
 
 return coords

def GetInteractionEnergy(coords,atoms_per_mol):
	
 cluster_size = len(coords)
 if ( cluster_size%atoms_per_mol !=0 ):
	print "Corrupt data"
	exit(2)
 cluster_size = int(cluster_size/atoms_per_mol)

 center = [0.0,0.0,0.0]
 energyELS = 0.0
 energyVDW = 0.0

 for imol in range(0,cluster_size):
	for jmol in range(0,cluster_size):
	 if (imol==jmol):
		continue
	 for iatom in range(0,atoms_per_mol):
		for jatom in range(0,atoms_per_mol):
		 rij = 0.0
		 for icoord in range(0,3):
			rij += ( coords[atoms_per_mol*imol+iatom][icoord] -
							 coords[atoms_per_mol*jmol+jatom][icoord] )**2
		 rij = sqrt(rij)
		 # electrostatic component
		 energyELS += charge[iatom] * charge[jatom] / (rij * ang2bohr)
		 # LJ component
		 if (iatom==0 and jatom==0):
			energyVDW += 4.0*epsilon*kjmol2hartree*( (sigma/rij)**12 - (sigma/rij)**6 )

 energy = [energyELS,energyVDW]

 return energy
	
# =============================================================
# ------------- the main program ------------------------------
# =============================================================
def main():

 # range of parameters (grid)
 # start value, end value, number of points
 ploops = [ 
						[2.0,8.0,31], # distance between oxygen atoms, do not set closer than 2 
						#[0.71*pi,2*pi,1], # 1st (xy-plane) mol rotation angle
						#[5.0*pi/4.0,2*pi,1] # 2nd (perp) mol rotation angle
						[0.0,pi,31], # 1st (xy-plane) mol rotation angle
						[0.0,2*pi,61] # 2nd (perp) mol rotation angle
					]

 atoms_per_mol = 3

 # end data section ====================================

 # deltas
 dp = [
				(ploops[0][1]-ploops[0][0])/(1 if ploops[0][2]==1 else ploops[0][2]-1),
				(ploops[1][1]-ploops[1][0])/(1 if ploops[1][2]==1 else ploops[1][2]-1),
				(ploops[2][1]-ploops[2][0])/(1 if ploops[2][2]==1 else ploops[2][2]-1)
			]

 params = [0.,0.,0.]
 toreturn = []

 for ip0 in range(0,ploops[0][2]):
	params[0] = ploops[0][0] + dp[0]*ip0 
	for ip1 in range(0,ploops[1][2]):
	 params[1] = ploops[1][0] + dp[1]*ip1
	 for ip2 in range(0,ploops[2][2]):
		params[2] = ploops[2][0] + dp[2]*ip2

		coords = GenerateCoords(params)
		energyCompon = GetInteractionEnergy(coords,atoms_per_mol)
		energy = energyCompon[0]+energyCompon[1]
		#print "%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f" % (params[0],params[1],params[2],energyCompon[0],energyCompon[1],energy)
		point = [params[0],params[1],params[2],energy]
		toreturn.append(point)

 print "Total points on a grid: %10d" % len(toreturn)

 return toreturn



