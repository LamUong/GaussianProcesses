from sklearn import datasets
from sklearn.gaussian_process import GaussianProcess
from sklearn.cross_validation import cross_val_score, KFold
import math
import numpy
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import random
import sys
import os
sys.path.append("/GaussianProcesses/tip3p/")
import datagenerator

def load_data(Quantum_enrgy,Parameters,trainingsize): 

	print("Loading data and generating descriptors")
	Data_array = []
	Energy_Array = []
	number = 1
	with open("datafile.x", "r") as ins:
		for line in ins:
			if(number%8==0):
				number = 0
			if number ==1:
				NewSubEnergyArray = [] 
			if number%8 != 0:
				if number ==7:
					Energy_Array.append(float(str(line).split()[0]))
					Data_array.append(NewSubEnergyArray)
				else:
					Splitted_Line = str(line).split()
					Coordinates = [float(Splitted_Line[1]),float(Splitted_Line[2]),float(Splitted_Line[3])]
					NewSubEnergyArray.append(Coordinates)
			number+=1

	New_energy_array_from_nonquantum_cal = []
	for coors in Data_array:
		energy = sum(datagenerator.GetInteractionEnergy(coors))
		New_energy_array_from_nonquantum_cal.append(energy)

	if Quantum_enrgy == True and Parameters == 'All_distances':
		x_train, y_train, x_test, y_test = get_parameters_all_distances(Data_array, Energy_Array,trainingsize)
	elif Quantum_enrgy == False and Parameters == 'All_distances':
		x_train, y_train, x_test, y_test = get_parameters_all_distances(Data_array, New_energy_array_from_nonquantum_cal,trainingsize)
	elif Quantum_enrgy == True and Parameters == 'Distances_and_Angles':
		x_train, y_train, x_test, y_test = get_old_parameters_Angles_and_Distances(Data_array, Energy_Array,trainingsize)
	elif Quantum_enrgy == False and Parameters == 'Distances_and_Angles':
		x_train, y_train, x_test, y_test = get_old_parameters_Angles_and_Distances(Data_array, New_energy_array_from_nonquantum_cal,trainingsize)
	return x_train, y_train, x_test, y_test


#Calculate Plane Equation given 3 coordinates
def plane_Equation(Coor1,Coor2,Coor3):
	v1 = (Coor3[0]-Coor1[0],Coor3[1]-Coor1[1],Coor3[2]-Coor1[2])
	v2 = (Coor2[0]-Coor1[0],Coor2[1]-Coor1[1],Coor2[2]-Coor1[2])
	crossv1v2 = (v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0])
	z = crossv1v2[0]*Coor1[0]+crossv1v2[1]*Coor1[1]+crossv1v2[2]*Coor1[2]
	return(crossv1v2[0],crossv1v2[1],crossv1v2[2],z)

#Calculate Dihedral Angle given 4 coordinates
def dihedral_Angle(Coor1,Coor2,Coor3,Coor4):
	Equ1 = plane_Equation(Coor1,Coor2,Coor3)
	Equ2 = plane_Equation(Coor2,Coor3,Coor4)
	a1 , a2 = Equ1[0],Equ2[0]
	b1 , b2 = Equ1[1],Equ2[1]
	c1 , c2 = Equ1[2],Equ2[2]
	d1 , d2 = Equ1[3],Equ2[3]
	cos = (a1*a2+b1*b2+c1*c2)/(math.sqrt(a1**2+b1**2+c1**2)*math.sqrt(a2**2+b2**2+c2**2))
	angle = math.acos(cos)
	angleDegree = math.degrees(angle)
	return angleDegree

#Caculating Cosine of angle at molA point
def cosine_from_coordinates(bisector, molA, molB):
	distanceA=float(distance_2_coordinates(bisector,molA))
	distanceB=float(distance_2_coordinates(molA,molB))
	distanceC=float(distance_2_coordinates(bisector,molB))
	cos = ((distanceB**2)+(distanceA**2)-(distanceC**2))/(2*(distanceA*distanceB))
	return cos

#Searching for the coordinate of Angle Bisector 
def recursive_search(Ratio,Coor2,Coor3,OriginCoor2,OriginCoo3):
	Mid = (float(Coor2[0]+Coor3[0])/2,float(Coor2[1]+Coor3[1])/2,float(Coor2[2]+Coor3[2])/2)
	D12 = distance_2_coordinates(OriginCoor2,Mid)
	D23 = distance_2_coordinates(Mid,OriginCoo3)
	CalculatedRatio = float(float(D12)/float(D23))
  
	if (CalculatedRatio-0.0001)<Ratio and (CalculatedRatio+0.00001)>Ratio:
		return Mid
	elif CalculatedRatio<Ratio:
		return recursive_search(Ratio,Mid,Coor3,OriginCoor2,OriginCoo3)
	elif CalculatedRatio>Ratio:
		return recursive_search(Ratio,Coor2,Mid,OriginCoor2,OriginCoo3)

#Searching for the coordinate of Angle Bisector 
def angle_bisector_coordinates(Coor1,Coor2,Coor3):
	#Coor1 is the apex
	D12 = distance_2_coordinates(Coor1,Coor2)
	D13 = distance_2_coordinates(Coor1,Coor3)
	Ratio = float(D12/D13)
	CoorBisect = recursive_search(Ratio,Coor2,Coor3,Coor2,Coor3)
	return CoorBisect

def distance_2_coordinates(Coor1,Coor2):
	return math.sqrt((Coor1[0]-Coor2[0])**2+(Coor1[1]-Coor2[1])**2+(Coor1[2]-Coor2[2])**2)

def average_between_distances(Listofdistances):
	Numpy_Array = numpy.asarray(Listofdistances) 
	return (numpy.sum(Numpy_Array)/len(Numpy_Array))

def symmetric_differences_of_2_distances(Distanceij,Distanceik):
	return ((Distanceij-Distanceik)**2)

def get_parameters_all_distances(data,energy,trainingsize):
	x_train = []
	y_train = []
	x_test = []
	y_test = []
	number = 0
	Listpairinter = [0,1,2,3,4,5]
	for dimer, ene in zip(data, energy):
		number+=1
		parameters = []
		for i in range(0,len(Listpairinter)-1,1):
			j=i+1
			while j <len(Listpairinter):
				parameters.append(distance_2_coordinates(dimer[i],dimer[j]))
				j+=1
		if number <trainingsize:
			x_test.append(parameters) 
			y_test.append(ene)
		else:
			x_train.append(parameters) 
			y_train.append(ene)
	return x_train, y_train, x_test, y_test 



def get_old_parameters_Angles_and_Distances(data,energy,trainingsize):
	x_train = []
	y_train = []
	x_test = []
	y_test = []
	number = 0
	
	for dimer, ene in zip(data, energy):
		number+= 1
		distance12 = distance_2_coordinates(dimer[0],dimer[1])
		distance13 = distance_2_coordinates(dimer[0],dimer[2])
		if distance12 < distance13:
			Hm = dimer[1]
		else:
			Hm = dimer[2]

		average1213 = average_between_distances([distance12,distance13])
		dif1213 = symmetric_differences_of_2_distances(distance12,distance13)
		distance23 = distance_2_coordinates(dimer[1],dimer[2])

		distance45 = distance_2_coordinates(dimer[3],dimer[4])
		distance46 = distance_2_coordinates(dimer[3],dimer[5])

		if distance45 < distance46:
			Hn = dimer[4]
		else:
			Hn = dimer[5]

		average4546 = average_between_distances([distance45,distance46])
		dif4546 = symmetric_differences_of_2_distances(distance45,distance46)
		distance56 = distance_2_coordinates(dimer[4],dimer[5])

		distance2Oxygens = distance_2_coordinates(dimer[0],dimer[3])
		Xa = angle_bisector_coordinates(dimer[0],dimer[1],dimer[2])
		Xb = angle_bisector_coordinates(dimer[3],dimer[4],dimer[5])
		cosXaO1O4 = cosine_from_coordinates(Xa, dimer[0], dimer[3])
		cosXbO4O1 = cosine_from_coordinates(Xb, dimer[3], dimer[0])
		diheXbO4O1Xa = dihedral_Angle(Xb,dimer[3],dimer[0],Xa)
		DiheHmXaO1O4 = dihedral_Angle(Hm,Xa,dimer[0],dimer[3])
		DiheHnXbO4O1 = dihedral_Angle(Hn,Xb,dimer[3],dimer[0])
	
		if number <trainingsize:
			x_test.append([average1213, dif1213, distance23, average4546, dif4546, distance56, 
				distance2Oxygens, cosXaO1O4, cosXbO4O1, diheXbO4O1Xa, DiheHmXaO1O4, DiheHnXbO4O1]) 
			y_test.append(ene)
		else:
			x_train.append([average1213, dif1213, distance23, average4546, dif4546, distance56, 
					distance2Oxygens, cosXaO1O4, cosXbO4O1, diheXbO4O1Xa, DiheHmXaO1O4, DiheHnXbO4O1]) 
			y_train.append(ene)
	
	return x_train, y_train, x_test, y_test 



def learning(regression,correlation,trainingsize,Quantum_enrgy,Parameters):
	x_train, y_train, x_test, y_test = load_data(Quantum_enrgy,Parameters,trainingsize)
	print len(x_train)
	print len(y_train)

	X_train, y_train = numpy.asarray(x_train) , numpy.asarray(y_train)
	X_val, y_val = numpy.asarray(x_test) , numpy.asarray(y_test)

	print("Generating GP")
	gp = GaussianProcess(corr=correlation, normalize=True, regr=regression, thetaL=1e-2, thetaU=1.0)
	gp.fit(X_train, y_train)
	return gp, x_train, y_train, x_test, y_test

def predict(index, gp, x_test, y_test):
	f2, MSE = gp.predict([x_test[index]], eval_MSE=True)
	predicted_energy = f2[0]
	validation_energy = y_test[index]
	sigma = np.sqrt(MSE)
	absolute_differences = math.fabs(predicted_energy-validation_energy)
	return sigma, predicted_energy, validation_energy, absolute_differences

def plot(x,y,e,truey):
	x = np.array(x)
	y = np.array(y) 
	e = np.array(e)
	truey = np.array(truey)
	f, (a0, a1) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[3, 1]})
	a0.plot(x, truey, 'ro')
	a0.errorbar(x, y, e, linestyle='None', marker='^')

	DayOfWeekOfCall = [1,2,3,4,5,6]
	DispatchesOnThisWeekday = [220, 220, 220, 220, 220, 220]

	LABELS = ["<3", "3-4", "4-5","5-6", "6-7", ">7"]
	a1.bar(y,x)

	a1.bar(DayOfWeekOfCall, DispatchesOnThisWeekday, align='center')
	a1.set_xticks(DayOfWeekOfCall)
	a1.set_xticklabels(LABELS)
	plt.savefig('pics/rzk.png')

def main(Quantum_enrgy, Parameters):
	'''
	regr = ['constant', 'linear', 'quadratic']
	corr = ['absolute_exponential', 'squared_exponential','generalized_exponential', 'cubic', 'linear']
	startlist = [1000,2000,3000,4000,5000,5200,5400,5600]
	nugget = 2.2204460492503131e-15
	'''
	regression = 'quadratic'
	correlation = 'squared_exponential'
	trainingsize = 200

	gp, x_train, y_train, x_test, y_test = learning(regression,correlation,trainingsize, Quantum_enrgy, Parameters)

	x = []
	y = []
	e = []
	truey = []
	print("Calculating MSE")
	if Quantum_enrgy == True and Parameters == 'All_distances':
		name = "Quantum,All_distances.txt"
	elif Quantum_enrgy == False and Parameters == 'All_distances':
		name = "Non_Quantum,All_distances.txt"
	elif Quantum_enrgy == True and Parameters == 'Distances_and_Angles':
		name = "Quantum,Distances_and_Angles.txt"
	elif Quantum_enrgy == False and Parameters == 'Distances_and_Angles':
		name = "Non_Quantum,Distances_and_Angles.txt"
	myfile = open(name, 'w')

	MSE =0
	MAPE = 0
	number =0
	for index in range(0,199,1):
		number+=1 
		print index
		sigma, predicted_energy, validation_energy, absolute_differences = predict(index, gp, x_test, y_test)
		MAPE += math.fabs(absolute_differences/validation_energy)
		MSE += absolute_differences**2
		
		x.append(x_test[index][6])
		y.append(predicted_energy)
		e.append(2*sigma)
		truey.append(validation_energy)
		myfile.write("%15.7f%15.7f%15.7f\n" % (predicted_energy,validation_energy,2*sigma) )

	MSE = MSE/trainingsize
	MAPE = MAPE*100/trainingsize
	print MSE
	print MAPE
	myfile.write("MSE = " + str(MSE))
	myfile.write("MAPE = " + str(MAPE))
	myfile.close()


	#plot(x,y,e,truey)

main(Quantum_enrgy = False, Parameters = 'Distances_and_Angles')