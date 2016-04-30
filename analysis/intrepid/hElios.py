#!/usr/bin/python
import numpy as np; import math; from scipy.constants import pi;

# set constant Radius
R=1
# number of iterations sampled from Gaussian to get statistics
numberOsamples=512
# total number of steps in the helix; NOTE: numstep+1 = number of points (atoms)
numstep=( '1/8', '1/4', '1/2', '1', '2', '4', '8')
# steps per turn of the helix; NOTE: stpturn/numsteps = number of turns
stpturn=( '1/256', '1/128', '1/64', '1/32', '1/16', '1/8', '1/4', '1/2', '1', '2', '4', '8', '16', '32', '64', '128', '256' )
# ratio of the Radius to the rise
Rrratio=( '1/256', '1/128', '1/64', '1/32', '1/16', '1/8', '1/4', '1/2', '1', '2', '4', '8', '16', '32', '64', '128', '256' )
# degree of std dev; values are percent of the mean value; a value of 12.5 would give a std dev that-s 12.5% of the mean
pctVari=( '1','3.125','6.25','12.5','25','50' )

# Do each sample-draw from Gaussian         (512)
for samples in range(numberOsamples):

  # Do total NUMber of STEPs in the helix   (7)
  for num in numstep:
  
    # Do STePs per TURN in the helix        (17)
    for stp in stpturn:
    
      # Do Radius-Rise Ratios               (17)
      for rrr in Rrratio:
  
        # Finally, we can start making helices
        ######################################
        
        # Make perfectly regular helix, so we can break it precisely below
        numturns = float(stpturn) / float(numstep)
        numatoms = int(numstep) + 1
        totsweep = numturns * 2 * pi
        steprise = R / rrr              # Radius / ratio_(Radius/rise) = rise
        heradius = R * rrr              # Radius * ratio_(Radius/rise) = Radius
        
        # First, let-s do PITCH
        for varians in pctVari:
          pitch = noise( varians, pitch ) 
          axis=0; radius=heradius; sweep=totsweep; nturns=numturns
          helix(axis, pitch, radius, sweep,  nturns )
          
        # Second, let-s do RADIUS
        for varians in pctVari:
          radius = noise( varians, radius )
          axis=0; pitch=pitch; sweep=totsweep; nturns=numturns
          helix(axis, pitch, radius, sweep, nturns )
      
        
