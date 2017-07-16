import sys, os, subprocess
from sys import argv
from math import exp, pi, cos, sin, sqrt, floor, ceil

#========= constants ============== 
ang2bohr = 1.889725989
hartree2kjmol=2625.49962
hartree2kcalmol=627.509469
kjmol2hartree=1.0/hartree2kjmol
kcalmol2hartree=1.0/hartree2kcalmol

atoms_per_mol = 3

radiusH = 0.4 * ang2bohr
radiusO = 0.5 * ang2bohr

# ==== all parameters must be in a.u. ====
# == TIP3P model of water
rOH_0 = 0.9572 * ang2bohr
aHOH_0 = 104.52 * pi / 180.0
charge = [-0.834, +0.417, +0.417] 
sigma = 3.15061 * ang2bohr
epsilon = 0.6364 * kjmol2hartree 
De_morse = 0.0
alpha = 0.0
k_r_aHOH = 0.0
k_rr = 0.0
k_aHOH = 0.0

# == flexible SCP model of water: JPC A, 108, 11056 (2004)
rOH_0 = 1.0 * ang2bohr
aHOH_0 = 109.47 * pi / 180.0
charge = [-0.82, +0.41, +0.41]
sigma = 3.166 * ang2bohr
epsilon = 0.1554 * kcalmol2hartree
De_morse = 101.9188 * kcalmol2hartree
alpha = 2.567 / ang2bohr 
k_r_aHOH = -211.4672 * kcalmol2hartree / ang2bohr**2
k_rr = 111.70765 * kcalmol2hartree / ang2bohr**2 
k_aHOH = 328.645606 * kcalmol2hartree / ang2bohr**2

def GenerateCoords(params,paramsID=0):
 
    coords = []
    
    nparams = len(params)
    for iparam in range(0,5):
        if (iparam == 0):
            # O-O distance
            R = params[0]
        elif (iparam == 1):
            # orientation of the molecule that lies in the H-bond plane
            # to get "standard" orientation, in which O..H-O angle is zero
            # Th must be equal to Pi - aHOH/2
            Th = params[1]
        elif (iparam == 2):
            # orientation of the molecule that does not lie in the H-bond plane
            Fi = params[2]
        elif (iparam == 3):
            if (paramsID==0):
                # covalent OH bond length
                if ( iparam < nparams ):
                    rOH = params[3]
                else:
                    rOH = rOH_0
            else:
                rOH = rOH_0
                # x coordinate of the proton
                if ( iparam < nparams ):
                    xHplus = params[3]
                else:
                    print ("x coord of the proton is missing")
                    exit(3)
        elif (iparam == 4):
            if (paramsID==0):
                # HOH intramolecular angle
                if ( iparam < nparams ):
                    aHOH = params[4]
                else:
                    aHOH = aHOH_0
            else:
                aHOH = aHOH_0
                # y coordinate of the proton
                if ( iparam < nparams ):
                    yHplus = params[4]
                else:
                    print ("y coord of the proton is missing")
                    exit(3)
    
    rOHxy = rOH * cos( aHOH/2.0 )
    rOHz  = rOH * sin( aHOH/2.0 )
    
    # generate 6 atoms
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
            if (paramsID==0):
             # use anlge
             ang = Th + aHOH/2.0
             x = R + rOH * cos( ang )
             y = rOH * sin( ang )
            else:
             # use proton cart coordinates
             x = xHplus
             y = yHplus
            xyz = [x,y,0.]
    
           else:
            print ("Atoms number is too large")
            exit(3)
    
           coords.append(xyz[:])
    
    return coords

def GetRadius(iatom):
 global atoms_per_mol
 if ( iatom%atoms_per_mol !=0 ):
     radius = radiusH
 else:
     radius = radiusO
 return radius

def CoordsAreReasonable(coords):

 natoms = len(coords)
 for iatom in range(0, natoms):
     radius1 = GetRadius(iatom)
     for jatom in range(iatom+1, natoms):
         radius2 = GetRadius(jatom)
         distance2 = 0.0
         for icoord in range(0, 3):
             distance2 += ( coords[iatom][icoord] - coords[jatom][icoord] )**2
         if ( distance2 < (radius1+radius2)**2 ):
             return False

 return True

# ==== coords must be in a.u. ====
def GetInteractionEnergy(coords):
	
 global atoms_per_mol
 cluster_size = len(coords)

 if ( cluster_size%atoms_per_mol !=0 ):
  print ("Corrupt data")
  exit(2)
 cluster_size = int(cluster_size/atoms_per_mol)

 energyELS = 0.0
 energyVDW = 0.0
 energyINTRA = 0.0

 for imol in range(0,cluster_size):

  # intramolecular terms
  rHH  = 0.0
  rOH1 = 0.0
  rOH2 = 0.0
  for icoord in range(0,3):
   rHH  += ( coords[atoms_per_mol*imol+1][icoord] - coords[atoms_per_mol*imol+2][icoord] )**2
   rOH1 += ( coords[atoms_per_mol*imol+0][icoord] - coords[atoms_per_mol*imol+1][icoord] )**2
   rOH2 += ( coords[atoms_per_mol*imol+0][icoord] - coords[atoms_per_mol*imol+2][icoord] )**2
  rHH = sqrt(rHH)
  rOH1 = sqrt(rOH1)
  rOH2 = sqrt(rOH2)
  
  rHH_0 = 2.0 * rOH1 * sin(aHOH_0/2.0) 
  
  delta_rHH  = rHH  - rHH_0
  delta_rOH1 = rOH1 - rOH_0
  delta_rOH2 = rOH2 - rOH_0
  
  energyINTRA += De_morse * ( (1.0-exp(alpha*delta_rOH1))**2 + (1.0-exp(alpha*delta_rOH2))**2 ) + k_r_aHOH * delta_rHH * (delta_rOH1 + delta_rOH2) + k_rr * delta_rOH1 * delta_rOH2 + 0.5 * k_aHOH * delta_rHH**2
  
  # intermolecular terms
  for jmol in range(imol,cluster_size):
   if (imol==jmol):
    continue
  	
   for iatom in range(0,atoms_per_mol):
    for jatom in range(0,atoms_per_mol):
  
     rij = 0.0
     for icoord in range(0,3):
       rij += ( coords[atoms_per_mol*imol+iatom][icoord] - coords[atoms_per_mol*jmol+jatom][icoord] )**2
     rij = sqrt(rij)
     # electrostatic component
     energyELS += charge[iatom] * charge[jatom] / rij
     # LJ component
     if (iatom==0 and jatom==0):
      energyVDW += 4.0*epsilon*( (sigma/rij)**12 - (sigma/rij)**6 )
  
 energy = [energyELS,energyVDW,energyINTRA]

 return energy
	
# =============================================================
# ------------- the main program ------------------------------
# =============================================================
def main(energyRequired):

 #if not energyRequired:
 #        ang2bohr = 1.0

 # range of parameters (grid)
 # start value, end value, number of points
 pgrid = [ 
		[3.0*ang2bohr,4.0*ang2bohr,3], # distance between oxygen atoms, do not set closer than 2 angstroms 
		#[0.66*pi,0.72*pi,3], # 1st (xy-plane) mol rotation angle
		#[5.0*pi/4.0,2*pi,1], # 2nd (perp) mol rotation angle
		#[0.0,pi,31], # 1st (xy-plane) mol rotation angle
		#[0.0,2*pi,61], # 2nd (perp) mol rotation angle
		[pi - aHOH_0/2.0, pi - aHOH_0/2.0, 1], # 1st (xy-plane) mol rotation angle
		[2.0*pi/3.0, 2.0*pi/3.0, 1], # 2nd (perp) mol rotation angle
		#[rOH_0*ang2bohr,rOH_0*ang2bohr,1] # INTRA-molecular OH distance
		#[aHOH_0,2*pi,1] # INTRA-molecular angle 
		[-2.*ang2bohr,6.*ang2bohr,81], # x coordinate of proton (along O-O direction)
		[-2.*ang2bohr,2.*ang2bohr,41] # y coordinate of proton (perp to O-O direction)
	  ]

 # end data section ====================================

 nparams = len(pgrid)
 params = [0.] * nparams
 datapoints = []

 GenerateGrid(0, pgrid, params, datapoints, energyRequired)

 print ("Number of parameters: %10d" % nparams)
 print ("Total points on a grid: %10d" % len(datapoints))

 return datapoints

def GenerateGrid(currentDim, pgrid, params, datapoints, energyRequired):

 ndims = len(pgrid)

 if ( currentDim < ndims):
  # delta
  dp = (pgrid[currentDim][1]-pgrid[currentDim][0])/(1 if pgrid[currentDim][2]==1 else pgrid[currentDim][2]-1)

  for idim in range(0,pgrid[currentDim][2]):
    params[currentDim] = pgrid[currentDim][0] + dp*idim
    GenerateGrid(currentDim+1, pgrid, params, datapoints, energyRequired)
 else:
    if (energyRequired):
        paramsID = 0
    else:
        paramsID = 1
    coords = GenerateCoords(params,paramsID)
    acceptCoords = CoordsAreReasonable(coords)
    if acceptCoords:
     if energyRequired:
      energyCompon = GetInteractionEnergy(coords)
      energy = sum(energyCompon)
      point = [params[:],energy]
     else:
      # energy will be evaluated using external (QC) package
      # based on the returned coordinates
      point = [params[:],coords[:]]
      #print point
      #print "%20.10f%20.10f%20.10f%20.10f" %(point[0][0], energyCompon[0],energyCompon[1],energyCompon[2])
      #print "%20.10f%20.10f" % (point[0][0], energy)
     datapoints.append(point[:])

