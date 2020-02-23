#!/usr/bin/env py -3

# Make a finite markov chain with transition probabilities
# Moving left (towards lower index) emits a 0, moving right emits a 1.
# Run the model to generate random bits.
#
# Maintains an occupancy count for each state.
#
# Occupancy distribution lets you compute the average entropy
# from the entropy of the bit emitted from each bit multiplied
# by the weight of the bin, determined by the occupancy.
#

import math  
import rdrand
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

r = rdrand.RdRandom()

chain_len = 100
x_low_limit = -1.0
x_high_limit = 1.0

pp = PdfPages('sigmoid_markov.pdf')

bitcount = 0
linecount = 0

def logistic_sigmoid(x):
    return 1.0/(1.0+math.exp(-x))

def hyperbolic_tangent_sigmoid(x):
    return (math.tanh(x))

def arctan_sigmoid(x):
    return (math.atan(x))

def gudermannian_sigmoid(x):
    return (2.0*math.atan(math.tanh(x/2.0)))

def erf_sigmoid(x):
    return (math.erf(x))

def algebraic_sigmoid(x):
    return (x/math.sqrt(1.0+(x**2)))
    
def output_bit(bit):
        global byte
        global linecount
        global bitcount
        byte = byte + (bit << bitcount)
        bitcount += 1
        if bitcount == 8:
            #print("%02X" % byte, end = '')
            bitcount = 0
            byte = 0
            linecount += 1
            if linecount == 32:
                linecount = 0
                #print("")
                    

def make_plot(chain_len, f, start_state):

    chain1 = [f(x_low_limit+(x*stepsize))  for x in range(chain_len)] 
    
    # Scale chain1 to 0-1
    minimum = min(chain1)
    maximum = max(chain1)
    
    h = maximum-minimum
    
    chain2 = [(x-minimum)/h for x in chain1] # scale and shift to 0-1 range
    
    chain = [(x,0) for x in chain2]          # add occupancy counts
    
    # Initialise markov chain
    #   leftp     = probability of moving to the left
    #   occupancy = count of the times that state was visited
        
    state = start_state
    (leftp,occupancy) = chain[state]
    
    # Generate 1 Mibibyte or 8 Mibibits of data
    # by running the markov model
    # and keeping record of the occupancy counts
     
    for i in range(8*(2**20)):
        if leftp > r.random():
            # Move left (decreasing index)
            if state > 0:
                state = state - 1
                bit = 0
        else:
            if state < (chain_len-1):
                state = state + 1
                bit = 1
        
        (leftp,occupancy) = chain[state]
        chain[state] = (leftp,occupancy+1)
        output_bit(bit)
    
    plist = list()
    for (leftp,occupancy) in chain:
        plist.append(leftp)
    
    occ_list = list()
    for (leftp,occupancy) in chain:
        occ_list.append(float(occupancy)/float(2**20))

    return plist, occ_list


start_state = (chain_len-1)//2

byte = 0

width = x_high_limit-x_low_limit
#print("width = ",width)


stepsize = (width/(chain_len-1))

functions = [   logistic_sigmoid,
                hyperbolic_tangent_sigmoid,
                arctan_sigmoid,
                gudermannian_sigmoid,
                erf_sigmoid,
                algebraic_sigmoid]

titles    = [   "logistic_sigmoid",
                "hyperbolic_tangent_sigmoid",
                "arctan_sigmoid",
                "gudermannian_sigmoid",
                "erf_sigmoid",
                "algebraic_sigmoid"]
                                
for f,title in zip(functions,titles):
    print(title)
    
    plist, occ_list = make_plot(chain_len, f, start_state)
        
    fig, ax = plt.subplots(1,2)
    fig.suptitle(title)
    
    ax[0].plot(plist)
    ax[1].plot(occ_list)
    
    pp.savefig()

pp.close()

    
#plt.plot(plist)
#plt.xlabel('state')
#plt.show()
#
#plt.plot(occ_list)
#plt.ylabel('State')
#plt.show()    
##print(float(occupancy)/(2**20))

