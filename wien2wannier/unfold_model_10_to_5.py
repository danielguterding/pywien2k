import sys
import os
import operator
import numpy as np

class HoppingElement:
  def __init__(self, x, y, z, o1, o2, t):
    self.x = float(x)
    self.y = float(y)
    self.z = float(z)
    self.o1 = int(o1)
    self.o2 = int(o2)
    self.t = float(t)
    self.abst = abs(self.t)
    self.v = np.array([self.x, self.y, self.z], dtype=float)
    
def get_input_elements(infilename):
  input_elements = []
  infilehandle = open(infilename, 'r')
  lines = infilehandle.readlines()
  infilehandle.close()
  
  for i,l in enumerate(lines):
    if(l.find("#x y z orbital1 orbital2 t") > -1):
      break
  lines = lines[i+1:] #take only lines with entries into account
  
  for l in lines:
    splitline = l.strip().split()
    x, y, z, o1, o2, t = splitline
    e = HoppingElement(x, y, z, o1, o2, t)
    input_elements.append(e)
    
  return input_elements
  
def group_elements_by_vector(input_elements):
  #sort input elements by z, y, x
  input_elements.sort(key=operator.attrgetter('z'))
  input_elements.sort(key=operator.attrgetter('y'))
  input_elements.sort(key=operator.attrgetter('x'))
  
  #make an array which stores information whether element is already grouped
  Nelements = len(input_elements)
  element_grouped = np.zeros(Nelements, dtype=int)
  
  #loop over elements and group them
  threshold = 1e-3
  grouped_elements = []
  firstungroupedidx = 0
  while(element_grouped.sum() < Nelements):
    newgroup = []
    #add first ungrouped element
    newgroup.append(input_elements[firstungroupedidx])
    element_grouped[firstungroupedidx] = 1
    
    #add elements to group by comparing hooping vectors
    for i in range(firstungroupedidx+1,Nelements):
      #check whether element is in group if it is currently ungrouped
      if(0 == element_grouped[i]):
	#check whether ungrouped element has the same hopping vector as first element in group
	if(np.linalg.norm(newgroup[0].v - input_elements[i].v) < threshold):
	  newgroup.append(input_elements[i])
	  element_grouped[i] = 1
	#assume that elements are ordered according to vectors so that we can stop searching once the first vector that does not match has been found
	else:
	  firstungroupedidx = i
	  break
    #append new group to list of groups
    grouped_elements.append(newgroup)
    
  return grouped_elements
  
def get_unfolded_elements(grouped_elements):
  np.set_printoptions(precision=2, linewidth=200)
  #build space group unfolding matrix according to Milan Tomic
  unfoldmatrix = np.identity(10)
  unfoldmatrix[0:5,5:10] = np.identity(5)
  unfoldmatrix[0,5] = -1.0 #assumes xz,yz orbitals with indices 0,1,5,6
  unfoldmatrix[1,6] = -1.0
  unfoldmatrix[5:10,0:5] = -unfoldmatrix[0:5,5:10].T
  #print unfoldmatrix
  
  #build t_ij matrix in each group and unfold it
  threshold = 1e-8
  unfolded_elements = []
  for group in grouped_elements:
    tmat = np.zeros((10,10))
    for e in group:
      tmat[e.o1,e.o2] = e.t
    #print tmat[0:5,0:5]
    #print tmat[5:10,5:10]
    #print 
    #if(np.linalg.norm(group[0].v) < 1e-7):
      #print tmat
    tmat = unfoldmatrix.dot(tmat.dot(unfoldmatrix.T))
    #if(np.linalg.norm(group[0].v) < 1e-7):
      #print
      #print tmat[0:5,0:5]
      #print tmat[5:10,5:10]
      #print
      #print tmat
      #print
    
    #take upper left block and save elements larger than threshold
    for o1 in range(5):
      for o2 in range(5):
	t = tmat[o1,o2]
	if(abs(t)>threshold):
	  e = HoppingElement(group[0].x, group[0].y, group[0].z, o1, o2, t)
	  unfolded_elements.append(e)

  return unfolded_elements
  
def write_output_file(outfilename, elements):
  outfilehandle = open(outfilename, 'w')
  outfilehandle.write("#hopping parameters in terms of the direct cell\n#x y z orbital1 orbital2 t\n")
  for e in elements:
    outfilehandle.write("% f % f % f %i %i % 1.14f\n" % (e.x, e.y, e.z, e.o1, e.o2, e.t))
  outfilehandle.close()
    
def main():
  if(len(sys.argv) == 3):
    infilename = sys.argv[1]
    outfilename = sys.argv[2]
    
    #read input model
    input_elements = get_input_elements(infilename)
    print 'Input model successfully read.'
    
    #group elements by real space hopping vector
    print 'Started grouping hopping elements.'
    grouped_elements = group_elements_by_vector(input_elements)
    print 'Finished grouping hopping elements.'
    
    #unfold element matrix for each group
    print 'Started unfolding elements.'
    unfolded_elements = get_unfolded_elements(grouped_elements)
    print 'Finished unfolding elements.'
    
    #sort unfolded elements by magnitude
    unfolded_elements.sort(key=operator.attrgetter('abst'))
    unfolded_elements.reverse()
    unfolded_elements.sort(key=operator.attrgetter('o2'))
    unfolded_elements.sort(key=operator.attrgetter('o1'))
    
    #write output file
    write_output_file(outfilename, unfolded_elements)
    print 'Output successfully written.'
    print 'Program finished.'
    
  else:
    print 'Wrong number of input parameters. Please supply input and output file names.'
  
main()