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
def load_data(npoints,ntests): 
	print("Generating data")

	data = datagenerator.main()

	print("Constructing test/train sets")
	random.shuffle(data)
	x = []
	y = []
	# reformat data
	for point in range(0,npoints):
		para1, para2, para3, energy = data[point]
		x.append([para1, para2, para3])
		y.append(energy)

	x_test  = x[0:ntests]
	y_test = y[0:ntests]
	x_train = x[ntests:npoints]
	y_train = y[ntests:npoints]

	return x_train, y_train, x_test, y_test


def learning(regression,correlation,x_train,y_train):

	print("Learning")

	X_train, Y_train = numpy.asarray(x_train) , numpy.asarray(y_train)

	gp = GaussianProcess(corr=correlation, normalize=True, regr=regression, thetaL=1e-2, thetaU=1.0)
	gp.fit(X_train, Y_train)
	return gp

def predict_point(gp, x_test):
	f2, MSE = gp.predict([x_test], eval_MSE=True)
	return f2[0], MSE

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

	npoints=10000
	ntests=1000

	x_train, y_train, x_test, y_test = load_data(npoints,ntests)
	print "Points in the training set: %10d" % len(y_train)
	print "Points in the     test set: %10d" % len(y_test)

	'''
	regr = ['constant', 'linear', 'quadratic']
	corr = ['absolute_exponential', 'squared_exponential','generalized_exponential', 'cubic', 'linear']
	startlist = [1000,2000,3000,4000,5000,5200,5400,5600]
	nugget = 2.2204460492503131e-15
	'''
	regression = 'quadratic'
	correlation = 'absolute_exponential'

	gp = learning(regression,correlation,x_train, y_train)

	print("Making predictions")
	x = []
	y = []
	e = []
	truey = []
	
	RMSE = 0.0
	MaxE = 0.0
	for index in range(ntests):
		predicted_energy, sigma2 = predict_point(gp, x_test[index])
		validation_energy = y_test[index]

                # statistics on the accuracy of the prediction
	        Error = predicted_energy-validation_energy
		absError = math.fabs(Error)
		RMSE += Error**2
	        sigma = np.sqrt(sigma2)
	
		if ( MaxE < absError ):
			MaxE = absError

		x.append(x_test[index][0]) # 1D plot as a function of one parameter
		y.append(predicted_energy)
		e.append(2*sigma)
		truey.append(validation_energy)

		print "%15.7f%15.7f%15.7f" % (predicted_energy,validation_energy,2*sigma) 

	RMSE = np.sqrt(RMSE/ntests)
	print "RMSE: %20.10f, MAXE: %20.10f" % (RMSE,MaxE)
	plot(x,y,e,truey)

main()
