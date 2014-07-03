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
  
def get_unfolded_elements(input_elements):
  unfolded_elements = []
  hrot = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]], np.float)
  signs = [1, -1, 1, -1, 1]
  # dxy dyz dx2y2 dxz dz2
  sites = np.array([[0.75, 0.25, 0.0], [0.25, 0.75, 0.0]], dtype=float)
  nbands = len(signs)
  orbstosites = np.array([0]*nbands + [1]*nbands, dtype=int)
  
  for e in input_elements:
    si = e.o1/nbands
    sj = e.o2/nbands
    bi = e.o1%nbands
    bj = e.o2%nbands
    
    sign = 0
    if si == 0 and sj == 0:
        sign = 1
    elif si == 1 and sj == 1:
        sign = signs[bi]*signs[bj]
    elif si == 0 and sj == 1:
        sign = signs[bj]
    elif si == 1 and sj == 0:
        sign = signs[bi]
    else:
      print 'Error in unfolding sign determination.'
        
    tvec = hrot.dot(e.v + sites[orbstosites[e.o2]] - sites[orbstosites[e.o1]])
    x = tvec[0]
    y = tvec[1]
    z = tvec[2]
    thop = 0.5*e.t*sign
    
    enew = HoppingElement(x, y, z, bi, bj, thop)
    unfolded_elements.append(enew)
  
  return unfolded_elements
    
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
  
def add_equivalent_elements(elements):
  #sort all elements by orbitals
  elements.sort(key=operator.attrgetter('o2'))
  elements.sort(key=operator.attrgetter('o1'))
  #sort all elements by coordinate
  elements.sort(key=operator.attrgetter('z'))
  elements.sort(key=operator.attrgetter('y'))
  elements.sort(key=operator.attrgetter('x'))
  
  elements_new = []
  #add all elements with equivalent hopping vectors and equal orbital indices
  lastungroupedidx=0
  nelements = len(elements)
  groupedelements=0
  while(groupedelements < nelements):
    group = [elements[lastungroupedidx]]
    groupedelements+=1
    e1 = group[0]
    for i in range(lastungroupedidx+1, nelements):
      e = elements[i]
      if((np.linalg.norm(e.v-e1.v) < 1e-7) and (e1.o1 == e.o1) and (e1.o2 == e.o2)):
	group.append(e)
	groupedelements+=1
      else:
	lastungroupedidx = i
	break
    #add all hopping elements in the group
    summedt = sum([e.t for e in group])
    enew = HoppingElement(e1.x, e1.y, e1.z, e1.o1, e1.o2, summedt)
    elements_new.append(enew)

  return elements_new
  
def write_output_file(outfilename, elements):
  outfilehandle = open(outfilename, 'w')
  outfilehandle.write("#hopping parameters in terms of the direct cell\n#x y z orbital1 orbital2 t\n")
  for e in elements:
    outfilehandle.write("% f % f % f %i %i % 1.14f\n" % (e.x, e.y, e.z, e.o1, e.o2, e.t))
  outfilehandle.close()
    
def main():
  if(len(sys.argv) == 4):
    hamilinfilename = sys.argv[1] #name of the Wannier90 Hamiltonian output case_hr.dat
    outfilename = sys.argv[2] #name of the Hamiltonian output file name
    elementcutoff = float(sys.argv[3]) #magnitude cutoff for hopping elements
    
    #get input elements with cutoff applied from Wannier90 file
    input_elements = get_input_elements(hamilinfilename, elementcutoff)
    print 'Input model successfully read.'
    
    #unfold hopping elements
    print 'Started unfolding elements.'
    unfolded_elements = get_unfolded_elements(input_elements)
    unfolded_elements = add_equivalent_elements(unfolded_elements)
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
    print 'Wrong number of input arguments. Please supply Wannier90 case_hr.dat, WIEN2k case.outputd, an ouptut file name and the numerical magnitude cutoff for the hopping elements.'
  
main()