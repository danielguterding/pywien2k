#Copyright (c) 2014, Daniel Guterding <guterding@itp.uni-frankfurt.de>
import os
import sys
import numpy as np

  def __init__(self, filename):
    print 'Do not use this class if your system is orthogonal. Wien2k only inteprets the k-vectors in the primitive cell if the lattice is orthogonal.'
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
    self.lattice_reciprocal_conventional = np.array(((l1,l2,l3),(l4,l5,l6),(l7,l8,l9)), dtype=float).transpose()
    
    while(-1 == lines[i].find('BR2_REC')):
      i+=1
    i+=1
    
    l1, l2, l3 = lines[i].split()
    i+=1
    l4, l5, l6 = lines[i].split()
    i+=1
    l7, l8, l9 = lines[i].split()
    i+=1
    self.lattice_reciprocal_primitive = np.array(((l1,l2,l3),(l4,l5,l6),(l7,l8,l9)), dtype=float).transpose()
    
    while(-1 == lines[i].find('BR1_DIR')):
      i+=1
    i+=1
    
    l1, l2, l3 = lines[i].split()
    i+=1
    l4, l5, l6 = lines[i].split()
    i+=1
    l7, l8, l9 = lines[i].split()
    i+=1
    self.lattice_direct_conventional = np.array(((l1,l2,l3),(l4,l5,l6),(l7,l8,l9)), dtype=float).transpose()
    
    while(-1 == lines[i].find('BR2_DIR')):
      i+=1
    i+=1
    
    l1, l2, l3 = lines[i].split()
    i+=1
    l4, l5, l6 = lines[i].split()
    i+=1
    l7, l8, l9 = lines[i].split()
    i+=1
    self.lattice_direct_primitive = np.array(((l1,l2,l3),(l4,l5,l6),(l7,l8,l9)), dtype=float).transpose()
	
    self.transform_direct = np.linalg.inv(self.lattice_direct_primitive).dot(self.lattice_direct_conventional)
    self.transform_direct_inverse = np.linalg.inv(self.transform_direct)
    
    self.transform_reciprocal = np.linalg.inv(self.lattice_reciprocal_primitive).dot(self.lattice_reciprocal_conventional)
    self.transform_reciprocal_inverse = np.linalg.inv(self.transform_reciprocal)

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
