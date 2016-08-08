from sklearn import datasets
from sklearn.gaussian_process import GaussianProcess
from sklearn.cross_validation import cross_val_score, KFold
import math
import numpy
from sys import argv

def planeDquation(Coor1,Coor2,Coor3):
	v1 = (Coor3[0]-Coor1[0],Coor3[1]-Coor1[1],Coor3[2]-Coor1[2])
	v2 = (Coor2[0]-Coor1[0],Coor2[1]-Coor1[1],Coor2[2]-Coor1[2])
	crossv1v2 = (v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0])
	z = crossv1v2[0]*Coor1[0]+crossv1v2[1]*Coor1[1]+crossv1v2[2]*Coor1[2]
	return(crossv1v2[0],crossv1v2[1],crossv1v2[2],z)
	
def dihedralAngle(Coor1,Coor2,Coor3,Coor4):
	Equ1 = planeDquation(Coor1,Coor2,Coor3)
	Equ2 = planeDquation(Coor2,Coor3,Coor4)
	a1 , a2 = Equ1[0],Equ2[0]
	b1 , b2 = Equ1[1],Equ2[1]
	c1 , c2 = Equ1[2],Equ2[2]
	d1 , d2 = Equ1[3],Equ2[3]
	cos = (a1*a2+b1*b2+c1*c2)/(math.sqrt(a1**2+b1**2+c1**2)*math.sqrt(a2**2+b2**2+c2**2))
	angle = math.acos(cos)
	angleDegree = math.degrees(angle)
	return angleDegree

def getting_cos(bisector, molA, molB):
	distanceA=float(distance_2_coordinates(bisector,molA))
 	distanceB=float(distance_2_coordinates(molA,molB))
 	distanceC=float(distance_2_coordinates(bisector,molB))
 	cos = ((distanceB**2)+(distanceA**2)-(distanceC**2))/(2*(distanceA*distanceB))
 	return cos

def recursivesearch(Ratio,Coor2,Coor3,OriginCoor2,OriginCoo3):
 	Mid = (float(Coor2[0]+Coor3[0])/2,float(Coor2[1]+Coor3[1])/2,float(Coor2[2]+Coor3[2])/2)
 	D12 = distance_2_coordinates(OriginCoor2,Mid)
 	D23 = distance_2_coordinates(Mid,OriginCoo3)
 	CalculatedRatio = float(float(D12)/float(D23))
  
 	if (CalculatedRatio-0.0001)<Ratio and (CalculatedRatio+0.00001)>Ratio:
   		return Mid
 	elif CalculatedRatio<Ratio:
 		return recursivesearch(Ratio,Mid,Coor3,OriginCoor2,OriginCoo3)
 	elif CalculatedRatio>Ratio:
 		return recursivesearch(Ratio,Coor2,Mid,OriginCoor2,OriginCoo3)

def angle_bisector_coordinates(Coor1,Coor2,Coor3):
  	#Coor1 is the apex
 	D12 = distance_2_coordinates(Coor1,Coor2)
 	D13 = distance_2_coordinates(Coor1,Coor3)
 	Ratio = float(D12/D13)
 	CoorBisect = recursivesearch(Ratio,Coor2,Coor3,Coor2,Coor3)
 	return CoorBisect

def distance_2_coordinates(Coor1,Coor2):
	return math.sqrt((Coor1[0]-Coor2[0])**2+(Coor1[1]-Coor2[1])**2+(Coor1[2]-Coor2[2])**2)

def average_between_distances(Listofdistances):
	Numpy_Array = numpy.asarray(Listofdistances) 
	return (numpy.sum(Numpy_Array)/len(Numpy_Array))

def symmetric_differences_of2_distances(Distanceij,Distanceik):
	return ((Distanceij-Distanceik)**2)

def Get_Parameters_for_Gaussian(data):
	ToReturn = []
	for dimer in data:

		distance12 = distance_2_coordinates(dimer[0],dimer[1])
		distance13 = distance_2_coordinates(dimer[0],dimer[2])
		if distance12 < distance13:
			Hm = dimer[1]
		else:
			Hm = dimer[2]

		average1213 = average_between_distances([distance12,distance13])
		dif1213 = symmetric_differences_of2_distances(distance12,distance13)
		distance23 = distance_2_coordinates(dimer[1],dimer[2])

		distance45 = distance_2_coordinates(dimer[3],dimer[4])
		distance46 = distance_2_coordinates(dimer[3],dimer[5])

		if distance45 < distance46:
			Hn = dimer[4]
		else:
			Hn = dimer[5]

		average4546 = average_between_distances([distance45,distance46])
		dif4546 = symmetric_differences_of2_distances(distance45,distance46)
		distance56 = distance_2_coordinates(dimer[4],dimer[5])

		distance2Oxygens = distance_2_coordinates(dimer[0],dimer[3])
		Xa = angle_bisector_coordinates(dimer[0],dimer[1],dimer[2])
		Xb = angle_bisector_coordinates(dimer[3],dimer[4],dimer[5])
		cosXaO1O4 = getting_cos(Xa, dimer[0], dimer[3])
		cosXbO4O1 = getting_cos(Xb, dimer[3], dimer[0])
		diheXbO4O1Xa = dihedralAngle(Xb,dimer[3],dimer[0],Xa)
		DiheHmXaO1O4 = dihedralAngle(Hm,Xa,dimer[0],dimer[3])
		DiheHnXbO4O1 = dihedralAngle(Hn,Xb,dimer[3],dimer[0])
		ToReturn.append([average1213, dif1213, distance23, average4546, dif4546, distance56, distance2Oxygens, cosXaO1O4, cosXbO4O1, diheXbO4O1Xa, DiheHmXaO1O4, DiheHnXbO4O1]) 
	
	return ToReturn

def load_data(): 
	Data_array = []
	Energy_Array = []
	number = 1
	with open("datafile.x", "r") as ins:
		for line in ins:
			if(number%8==0):
				number = 0
			if number ==1:
				newsubarray = [] 
			if number%8 != 0:
				if number ==7:
					Energy_Array.append(float(str(line).split()[0]))
					Data_array.append(newsubarray)
				else:
					splitted_Line = str(line).split()
					Coordinates = (float(splitted_Line[1]),float(splitted_Line[2]),float(splitted_Line[3]))
					newsubarray.append(Coordinates)
			number+=1
	parameters = Get_Parameters_for_Gaussian(Data_array)
	return parameters,Energy_Array

def learning():
	print("Generating GP")
	Dataarray,Energyarray = load_data()
	X_train, y_train = numpy.asarray(Dataarray[1000:]) , numpy.asarray(Energyarray[1000:])
	gp = GaussianProcess(corr='squared_exponential', normalize=True, regr='linear', thetaL=1e-2, thetaU=1.0)
	gp.fit(X_train, y_train)
	return gp, Dataarray, Energyarray

def main(x, gp, Dataarray, Energyarray):
	f2, MSE2 = gp.predict([Dataarray[x]], eval_MSE=True)
	print f2
	print Energyarray[x]
	print MSE2

gp, Dataarray, Energyarray = learning()

while True:
	Index = raw_input("List Index? ")
	main(int(Index),gp,Dataarray,Energyarray)
	print("\n")
