#Copyright (c) 2014, Daniel Guterding <guterding@itp.uni-frankfurt.de>
import os
import sys
import numpy as np

class Wien2kConventionalToPrimitive:
  def __init__(self, filename):
    filehandle = open(filename, 'r')
    lines = filehandle.readlines()
    filehandle.close()
    
    self.lattice_direct_conventional = np.zeros((3,3), dtype=float)
    self.lattice_direct_primitive = np.zeros((3,3), dtype=float)
    self.lattice_reciprocal_conventional = np.zeros((3,3), dtype=float)
    self.lattice_reciprocal_primitive = np.zeros((3,3), dtype=float)
    
    i = 0
    while(-1 == lines[i].find('BR1_REC')):
      i+=1
    i+=1
    
    l1, l2, l3 = lines[i].split()
    i+=1
    l4, l5, l6 = lines[i].split()
    i+=1
    l7, l8, l9 = lines[i].split()
    i+=1
    self.lattice_reciprocal_conventional = np.array(((l1,l2,l3),(l4,l5,l6),(l7,l8,l9)), dtype=float)
    
    while(-1 == lines[i].find('BR2_REC')):
      i+=1
    i+=1
    
    l1, l2, l3 = lines[i].split()
    i+=1
    l4, l5, l6 = lines[i].split()
    i+=1
    l7, l8, l9 = lines[i].split()
    i+=1
    self.lattice_reciprocal_primitive = np.array(((l1,l2,l3),(l4,l5,l6),(l7,l8,l9)), dtype=float)
    
    while(-1 == lines[i].find('BR1_DIR')):
      i+=1
    i+=1
    
    l1, l2, l3 = lines[i].split()
    i+=1
    l4, l5, l6 = lines[i].split()
    i+=1
    l7, l8, l9 = lines[i].split()
    i+=1
    self.lattice_direct_conventional = np.array(((l1,l2,l3),(l4,l5,l6),(l7,l8,l9)), dtype=float)
    
    while(-1 == lines[i].find('BR2_DIR')):
      i+=1
    i+=1
    
    l1, l2, l3 = lines[i].split()
    i+=1
    l4, l5, l6 = lines[i].split()
    i+=1
    l7, l8, l9 = lines[i].split()
    i+=1
    self.lattice_direct_primitive = np.array(((l1,l2,l3),(l4,l5,l6),(l7,l8,l9)), dtype=float)
	
    self.transform_direct = self.lattice_direct_primitive.dot(np.linalg.inv(self.lattice_direct_conventional))
    self.transform_direct_inverse = self.lattice_direct_conventional.dot(np.linalg.inv(self.lattice_direct_primitive))
    
    self.transform_reciprocal = self.lattice_reciprocal_primitive.dot(np.linalg.inv(self.lattice_reciprocal_conventional))
    self.transform_reciprocal_inverse = self.lattice_reciprocal_conventional.dot(np.linalg.inv(self.lattice_reciprocal_primitive))
    
    #print self.lattice_reciprocal_conventional
    #print self.lattice_reciprocal_primitive
    #print self.lattice_direct_conventional
    #print self.lattice_direct_primitive
    
    #print np.linalg.inv(self.lattice_reciprocal_conventional)
    #print np.linalg.inv(self.lattice_reciprocal_primitive)
    #print np.linalg.inv(self.lattice_direct_conventional)
    #print np.linalg.inv(self.lattice_direct_primitive)
    
    #print self.transform_direct
    #print self.transform_direct_inverse
    #print self.transform_reciprocal
    #print self.transform_reciprocal_inverse

  def direct_primitive_to_conventional(self, vec):
    vec = np.array(vec, dtype=float)
    return self.transform_direct_inverse.dot(vec)
    
  def direct_conventional_to_primitive(self, vec):
    vec = np.array(vec, dtype=float)
    return self.transform_direct.dot(vec)
    
  def reciprocal_primitive_to_conventional(self, vec):
    vec = np.array(vec, dtype=float)
    return self.transform_reciprocal_inverse.dot(vec)
    
  def reciprocal_conventional_to_primitive(self, vec):
    vec = np.array(vec, dtype=float)
    return self.transform_reciprocal.dot(vec)
    
def main():
  
  if(1<len(sys.argv)):
    filename = str(sys.argv[1])
    mode = str(sys.argv[2])
    x = float(sys.argv[3])
    y = float(sys.argv[4])
    z = float(sys.argv[5])
  
    converter = Wien2kConventionalToPrimitive(filename)
    if('dc2p' == mode):
      print converter.direct_conventional_to_primitive([x,y,z])
    elif('dp2c' == mode):
      print converter.direct_primitive_to_conventional([x,y,z])
    elif('rc2p' == mode):
      print converter.reciprocal_conventional_to_primitive([x,y,z])
    elif('rp2c' == mode):
      print converter.reciprocal_primitive_to_conventional([x,y,z])
    else:
      print 'Warning! Mode not recognized.'
  else:
    print "No input file supplied! Please pass the path to case.outputd as command line argument."
  
main()
