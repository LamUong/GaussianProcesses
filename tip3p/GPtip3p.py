from sklearn import datasets
from sklearn.gaussian_process import GaussianProcess
from sklearn.cross_validation import cross_val_score, KFold
import math
import numpy
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import datagenerator
import random
'''
I only used 10 000 data for training caused it was taking a long time (several minutes). The result was plotted and printed. 
The accuracy was quite good at 10000 training. It can be further increased. MSE = 5E-7
 
I also tried 15 000 data point. The MSE still decreased by half compared to 10 000. The accuracy further increased. MSE = 2.5E-7

'''
def load_data(): 
	data = datagenerator.main()
	random.shuffle(data)
	x = []
	y = []
	for point in data:
		para1, para2, para3, energy = point
		x.append([para1, para2, para3])
		y.append(energy)

	x_train = x[200:15000]
	y_train = y[200:15000]
	x_test  = x[0:200]
	y_test = y[0:200]

	return x_train, y_train, x_test, y_test


def learning(regression,correlation):
	print("Loading data and generating descriptors")
	x_train, y_train, x_test, y_test = load_data()
	print len(x_train)
	print len(y_train)

	X_train, Y_train = numpy.asarray(x_train) , numpy.asarray(y_train)

	print("Generating GP")
	gp = GaussianProcess(corr=correlation, normalize=True, regr=regression, thetaL=1e-2, thetaU=1.0)
	gp.fit(X_train, Y_train)
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

	f, (a0) = plt.subplots(1,1)
	a0.plot(x, truey, 'ro')
	a0.errorbar(x, y, e, linestyle='None', marker='^')
	plt.savefig('rzk.png')

def main():
	'''
	regr = ['constant', 'linear', 'quadratic']
	corr = ['absolute_exponential', 'squared_exponential','generalized_exponential', 'cubic', 'linear']
	startlist = [1000,2000,3000,4000,5000,5200,5400,5600]
	nugget = 2.2204460492503131e-15
	'''
	regression = 'quadratic'
	correlation = 'absolute_exponential'

	gp, x_train, y_train, x_test, y_test = learning(regression,correlation)

	x = []
	y = []
	e = []
	truey = []
	print("Calculating MSE")
	
	MSE =0
	number =0
	for index in range(0,200,1):
		number+=1 
		print index
		sigma, predicted_energy, validation_energy, absolute_differences = predict(index, gp, x_test, y_test)

		MSE += absolute_differences**2
		
		x.append(x_test[index][0])
		y.append(predicted_energy)
		e.append(2*sigma)
		truey.append(validation_energy)
		print "%15.7f%15.7f%15.7f" % (predicted_energy,validation_energy,2*sigma) 

	MSE = MSE/200
	print MSE
	plot(x,y,e,truey)

main()
