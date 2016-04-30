#!/usr/bin/python
import numpy as np; import math; from scipy.constants import pi; from fractions import Fraction

# set constant Radius
R=1  #we test many "radius by scaling this using Rrratio...
pitch=1  #we need a Rkratio below, like Rrratio... now it-s just ONE pitch for code testing...

# number of iterations sampled from Gaussian to get statistics
numberOsamples=2

# steps per turn of the helix; NOTE: stpturn/numsteps = number of turns
stpturn=( '1/8', '1/4' ) #, '1/2', '1', '2', '4', '8' )

# ratio of the Radius to the rise
Rrratio=( '1/8', '1/4' ) #, '1/2', '1', '2', '4', '8' )

# degree of std dev; values are percent of the mean value; a value of 12.5 would give a std dev that-s 12.5% of the mean
pctVari=( '1', '3.125' )
#-----------------------------------------------------------------------------------
#  Define noise machine...
def noise(pctV, iparm):
  V=float(pctV)/float(100)
  mean_of_distro = float(iparm)
  std_dev_distro = float(V)*float(iparm) #stdev of distro is the percent of mean

  # CHECK to make sure i get what i want:
  #print 'V: %s' % V, 'mean_of_distro: ', mean_of_distro, 'std_dev_distro: ', std_dev_distro
 
  # Damien --- is the abs(std_dev_distro) below legal? it-s the variance that drives Gaussian (stdev**2), so i think this is ok...
  #  i did abs() becuz i was getting negative stdev-s, i think because the 
  oparm=np.random.normal(mean_of_distro,np.abs(std_dev_distro),size=None) 
  return oparm
#-----------------------------------------------------------------------------------
#  Define helix machine...
def helix(axis,pitch,radius,sweep,atoms):
  x=[];y=[];z=[]
  twist=sweep/atoms; t=0
  for pt in range(atoms):
    t=t+twist 
    x.append( radius * np.sin( np.radians(t) ) ); y.append( radius * np.cos( np.radians(t) ) ); z.append( pitch  * np.radians(t)
  return x,z,y #return side-ways helix b/c when phi=0 in spherical coordinates, things get weird.
#-----------------------------------------------------------------------------------
# Do total NUMber of STEPs in the helix   (7)
x=[]; y=[]; z=[]
for num in range(3,8,1):
  nu = int(num)
  # Do STePs per TURN in the helix        (17)
  for stp in stpturn:
    # Do Radius-Rise Ratios               (17)
    for rrr in Rrratio:
      # Finally, we can start making helices
      ######################################
      # Convert fractions in the lists to floats...
      if '/' in stp:
        st = float(Fraction(stp))
      else:
        st = float(stp)
      if '/' in rrr:
        rr = float(Fraction(rrr))
      else:
        rr = float(rrr)
      # Make perfectly regular helix, so we can break it precisely below
      numturns = float(st) / float(nu)
      numatoms = int(nu) + 1
      totsweep = numturns * 2 * pi
      steprise = R / rr              # Radius / ratio_(Radius/rise) = rise
      heradius = R * rr              # Radius * ratio_(Radius/rise) = Radius

      # CHECKto make sure what i get is what i want
      print 'num_atoms: %s' % numatoms, 'total_steps: ', num, 'steps_per_turn: ', stp, 'Radius-rise_ratio', rrr
      # First, let-s do PITCH
      for varians in pctVari:
        # Do the samples down here...
        out=''
        for samples in range(numberOsamples):
          pitch = noise( varians, pitch ) 
          axis=0; radius=heradius; sweep=totsweep; nturns=numatoms;
          x,z,y = helix(axis, pitch, radius, sweep,  nturns )
          out=out+str(x)+' '+str(z)+' '+str(y)+'\n'
#turn-off file IO for now...
"""
        with open('output_num'+str(nu)+'.stp'+str(st)+'.Rrr'+str(rr)+'.var'+str(varians), 'w') as f:
          f.write(out)
"""
          
        # Second, let-s do RADIUS
"""        for varians in pctVari:
          radius = noise( varians, radius )
          axis=0; pitch=pitch; sweep=totsweep; nturns=numatoms; #nturns=numturns
          helix(axis, pitch, radius, sweep, nturns )
"""
