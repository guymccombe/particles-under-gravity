import sys
from random import random


def createParticles( numberOfParticles, minMass, maxMass ):
  numberOfParticlesPerAxis = (int)(round(numberOfParticles**(1.0/3.0)))
  particleString = ""
  print( "create " + str(numberOfParticlesPerAxis) + " particles per axis" )
  # h = 1.0 /numberOfParticlesPerAxis
  h = 1.0
  for x in range(0,numberOfParticlesPerAxis):
   for y in range(0,numberOfParticlesPerAxis):
    for z in range(0,numberOfParticlesPerAxis):
      xPos = (random() - 0.5)*0.9*h + x * h
      yPos = (random() - 0.5)*0.9*h + y * h
      zPos = (random() - 0.5)*0.9*h + z * h
      mass = random()*(maxMass-minMass) + minMass
      particleString = particleString + " " + str(xPos) + " " + str(yPos) + " " + str(zPos) + " 0 0 0 " + str(mass) + "     "
  return particleString


def dumpToFile(particleString):       
  dumpFile = open( "initial-conditions.txt", "w" )
  dumpFile.write( particleString )
  dumpFile = open( "initial-conditions-with-plot.txt", "w" )
  dumpFile.write( "0.001  1.0  0.00001 " + particleString )
  dumpFile = open( "initial-conditions-no-plot.txt", "w" )
  dumpFile.write( "0.0  1.0  0.00001 " + particleString )


def informUser(particleString):
  if particleString.__len__()<80:
    print( particleString )
  print( "Written the following files:" )
  print( "  initial-conditions.txt\tPlain initial conditions")       
  print( "  initial-conditions-with-plot.txt\tPlain initial conditions")       
  print( "  initial-conditions-no-plot.txt\tPlain initial conditions")       
  print( "")       
  print( "cat initial-conditions.txt\t\tShow initial condition")
  print( "./mysimulator $(cat initial-conditions-with-plot.txt)")       



if __name__ =="__main__":
  if len(sys.argv)!=4:
    print( "Simple tool to create a random number of particles with random properties ")
    print( "within a cube. The particles are rougly spaced along a 1x1x1 voxel layout.")
    print( "")
    print( "Usage: python ./" + sys.argv[0] + " number-of-particles  min-mass  max-mass")
    print( "" )
    print( "Example: python create-initial-conditions.py 170 0.04 0.8" )
  else:
    numberOfParticles = (int)(sys.argv[1])
    minMass           = (float)(sys.argv[2])
    maxMass           = (float)(sys.argv[3])
    particleString = createParticles( numberOfParticles, minMass, maxMass )
    informUser(particleString)
    dumpToFile(particleString)
  
  