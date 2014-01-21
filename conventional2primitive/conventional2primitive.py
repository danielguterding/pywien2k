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
    while(-1 == lines[i].find('LATTICE CONSTANTS ARE')):
      i+=1
    a, b, c = lines[i].split()[-3:]
    
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
    
    
    #divide all matrices by lattice constants in order to go to reduced coordinates in the respective cell
    v = np.zeros(3, dtype=float)
    v = np.array((a,b,c), dtype=float)
    for i in range(3):
      for j in range(3):
	self.lattice_reciprocal_conventional[i][j] *= v[j]/2.0/np.pi
	self.lattice_reciprocal_primitive[i][j] *= v[j]/2.0/np.pi
	self.lattice_direct_conventional[i][j] /= v[j]
	self.lattice_direct_primitive[i][j] /= v[j]
    
    #print self.lattice_reciprocal_conventional
    #print self.lattice_reciprocal_primitive
    
    #print self.lattice_direct_conventional
    #print self.lattice_direct_primitive
    
    #print self.lattice_reciprocal_conventional * np.linalg.inv(self.lattice_reciprocal_primitive).transpose()
    #print self.lattice_direct_conventional * np.linalg.inv(self.lattice_direct_primitive).transpose()

  def direct_primitive_to_conventional(self, vec):
    vec = np.array(vec, dtype=float)
    return np.linalg.inv(self.lattice_direct_primitive).dot(vec)
    
  def direct_conventional_to_primitive(self, vec):
    vec = np.array(vec, dtype=float)
    return self.lattice_direct_primitive.dot(vec)
    
  def reciprocal_primitive_to_conventional(self, vec):
    vec = np.array(vec, dtype=float)
    return np.linalg.inv(self.lattice_reciprocal_primitive).dot(vec)
    
  def reciprocal_conventional_to_primitive(self, vec):
    vec = np.array(vec, dtype=float)
    return self.lattice_reciprocal_primitive.dot(vec)
    
def main():
  
  if(1<len(sys.argv)):
    filename = str(sys.argv[1])
  
    converter = Wien2kConventionalToPrimitive(filename)
    print converter.direct_conventional_to_primitive([1,0,0])
    print converter.direct_primitive_to_conventional([0,0,1])
    print converter.reciprocal_conventional_to_primitive([0,0,1])
    print converter.reciprocal_primitive_to_conventional([0,0,2])
  else:
    print "No input file supplied! Aborting. Please pass the path to case.outputd."
  
main()