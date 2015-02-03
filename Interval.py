from __future__ import division
from scipy.stats import binom
import random
import math
import numpy as np
import matplotlib.pyplot as plt
"""
In the exercise above, you tried to find an intuitive cutoff for likelihood ratio
scores that would give you a reasonable interval in which to find the true value of
p. Now, we will empirically determine one way to construct such an interval. To do
so, we will ask how far away from the true value of a parameter the ML estimate
might stray. Use this procedure: (1) start with a known value for p, (2) simulate
a bunch of datasets, (3) find ML parameter estimates for each simulation, and then
(4) calculate the likelihood ratios comparing the true parameter values and the ML
estimates. When you do this, you will be constructing a null distribution of
likelihood ratios that might be expected if the value of p you picked in (1)
was true. Note that the ML values for these replicates are very often greater than
L(true value of P), because the ML value can only ever be >= L(true value). Once
you have this distribution, find the likelihood ratio cutoff you need to ensure
that the probability of seeing an LR score that big or greater is <= 5%.
"""
def multiplier(num1, num2=1):
	""" Multiplies all consecutive numbers between num1 and num2 """
	start = max(num1, num2) #Get larger number
	stop = min(num1, num2) #Get smaller number
	product = stop #Initialize product as 1
	while int(start) != int(stop): #Loop while larger number is not equal to one less than smaller number
		product *= start #Multiply the results
		start -= 1 #Decrease larger number by 1
	return product #Return product
		
def bincoeffb(k, n):
	""" Calculate binomial coefficient (factoring first) from n sets choosing k subsets"""
	try:
		coeff = multiplier(n, n-k+1)/multiplier(k)
	except:
		coeff = 0
	return coeff
	
def bernoulli(k, n, p):
	"""Calculate a Probability Mass Function (PMF) of a binomial distribution with k successes in n Bernoulli trials with probability p"""
	PMF = bincoeffb(k, n)*pow(p,k)*pow((1-p),(n-k))
	return PMF

def hillclimb(data, prob, likeli=0, diff = 0.1, oldvalue=1, direction="RIGHT"):
	"""Calculate maximum likelihood for a dataset"""
	while diff >= 0.001: #Terminator loop
		if likeli == oldvalue: #Terminate recursion if values equal and change difference parameter to start new recursive list
			likeli = likeli
			prob = prob
			diff = diff/5
			oldvalue = 1
			hillclimb(data, prob=prob, likeli=likeli, diff=diff, oldvalue=oldvalue, direction="RIGHT")
		else: #Recurvise function to look for maximum point
			oldvalue = likeli
			problist = [prob-diff, prob, prob+diff]
			valuelist = [bernoulli(data[0], data[1], prob-diff), bernoulli(data[0], data[1], prob), bernoulli(data[0], data[1], prob+diff)]
			if direction == "LEFT":
				problist = problist[::-1]
				valuelist = valuelist[::-1]
			likeli = max(valuelist) #Choose maximum value of likelihood
			prob = problist[valuelist.index(likeli)]
			if prob<=0: #The program always goes left if two values are equal due to the max function. THis shifts it right.
				direction = "LEFT"
			hillclimb(data, prob=prob, likeli=likeli, diff=diff, oldvalue=oldvalue, direction=direction)
	return likeli, prob
	
# Set a starting, true value for p
trueP = 0.5

# Simulate 1,000 datasets of 200 trials from a binomial with this p
datalist = []
trials = 200
for i in range(1000):
	data = binom.rvs(trials,trueP)
	datalist.append(data)
	

# Now find ML parameter estimates for each of these trials
MLvals = []
likelivals = []
for val in datalist:
	likes = bernoulli(val, trials, trueP)
	maxlikeli, probability = hillclimb((val, trials), 0.5)
	likelivals.append(likes)
	MLvals.append(maxlikeli)

# Calculate likelihood ratios comparing L(trueP) in the numerator to the maximum likelihood (ML) in the denominator. Sort the results and find the value corresponding to the 95th percentile.
likeliratio = sorted(map(lambda x: likelivals[MLvals.index(x)]/x, MLvals))
print np.percentile(likeliratio, 95)

# Now, convert the likelihood ratios (LRs) to -2ln(LRs) values. Find the 95th percentile of these values. Compare these values to this table:
likelilog = map(lambda x: -2*math.log(x), likeliratio)
plt.hist(likelilog)
plt.show()
print np.percentile(likelilog, 95)

# https://people.richland.edu/james/lecture/m170/tbl-chi.html. In particular, look at the 0.05 column. Do any of these values seem similar to the one you calculated? Any idea why that particular cell would be meaningful?
# Based on your results (and the values in the table), what LR statistic value [-2ln(LR)] indicates that a null value of p is far enough away from the ML value that an LR of that size is <=5% probable if that value of p was true?
# Using this cutoff, what interval might you report for the 5- and 20-trial data sets above?
# We've talked in previous classes about two ways to interpret probabilities. Which interpretation are we using here to define these intervals?
