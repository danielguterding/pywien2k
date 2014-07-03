import sys
import os
import operator
import numpy as np

class HoppingElement:
  def __init__(self, x, y, z, o1, o2, t):
    self.x = int(x)
    self.y = int(y)
    self.z = int(z)
    self.o1 = int(o1)
    self.o2 = int(o2)
    self.t = float(t)
    self.abst = abs(self.t)
    self.v = np.array([self.x, self.y, self.z], dtype=float)
  def set_vector(self, v):
    self.v = v
    self.x = self.v[0]
    self.y = self.v[1]
    self.z = self.v[2]
    
def get_input_elements(infilename, elementcutoff):
  input_elements = []
  infilehandle_elements = open(infilename, 'r')
  for line in infilehandle_elements:
    splitline = line.strip().split()
    if(len(splitline) == 7 and not "1" == splitline[-1]):
      x, y, z, o1, o2, t, dt = splitline[0:7]
      o1 = int(o1)-1
      o2 = int(o2)-1
      entry = HoppingElement(x, y, z, o1, o2, t)
      if(entry.abst > elementcutoff):
	input_elements.append(entry)
  infilehandle_elements.close()
  
  return input_elements

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
  
def write_output_file(outfilename, elements):
  outfilehandle = open(outfilename, 'w')
  outfilehandle.write("#hopping parameters in terms of the direct cell\n#x y z orbital1 orbital2 t\n")
  for e in elements:
    outfilehandle.write("% f % f % f %i %i % 1.14f\n" % (e.x, e.y, e.z, e.o1, e.o2, e.t))
  outfilehandle.close()
  
def main():
  if(len(sys.argv) == 5):
    hamilinfilename = sys.argv[1] #name of the Wannier90 Hamiltonian output case_hr.dat
    outputdinfilename = sys.argv[2] #name of the Wien2k file case.outputd 
    outfilename = sys.argv[3] #name of the Hamiltonian output file name
    elementcutoff = float(sys.argv[4]) #magnitude cutoff for hopping elements
    
    #get input elements with cutoff applied from Wannier90 file
    input_elements = get_input_elements(hamilinfilename, elementcutoff)
    
    #get unit cell conversion matrices from WIEN2k file
    cell_converter = Wien2kConventionalToPrimitive(outputdinfilename)
    
    #convert all hopping vectors from primitive to conventional cell
    converted_elements = []
    for e in input_elements:
      e.set_vector(cell_converter.direct_primitive_to_conventional(e.v))
      converted_elements.append(e)
    input_elements = []
    
    #sort entries by magnitude
    converted_elements.sort(key=operator.attrgetter('abst'))
    converted_elements.reverse()
    converted_elements.sort(key=operator.attrgetter('o2'))
    converted_elements.sort(key=operator.attrgetter('o1'))
    
    #write output file
    write_output_file(outfilename, converted_elements)
    
    print 'Model successfully converted. %i hopping elements remaining.' % len(converted_elements)

  else:
    print 'Wrong number of input arguments. Please supply Wannier90 case_hr.dat, WIEN2k case.outputd, an ouptut file name and the numerical magnitude cutoff for the hopping elements.'
  
main()
