
############# SIMULATOR OF THE COHERENT RECEIVER CHANNELS ############
####         with random-walk functions for delta, theta, phi   #####  


from numpy import *
import numpy.random as npr


#####################################
###########----entries----###########



class RandomWalk:

	last=0
	
	def funrand(self,t):

		step = npr.normal(0,1)
		self.last+= step #+ (10**(-5))*npr.normal())	
		return self.last

######################################




