import os
import math
import numpy as np

class bandweights:
  def __init__(self, bandnumber):
    self.bn = bandnumber
    self.dtot = []
    self.dz2 = []
    self.dx2y2 = []
    self.dxy = []
    self.dxz = []
    self.dyz = []

def getqtldata(qtlfilename, orbitalmatching, atomnumber=1, interestingbands=[1]):
  filehandle = open(qtlfilename, 'r')
  lines = filehandle.readlines()
  filehandle.close()
  
  bandindices = []
  for i,l in enumerate(lines):
    if(l.find('BAND') > -1):
      bandindices.append([i+1])
      if(len(bandindices)>1):
	bandindices[-2].append(i-1)
  bandindices[-1].append(len(lines)-1)
  
  weightdata = []
  for i in interestingbands:
    weights = bandweights(i)
    for l in lines[bandindices[i][0]:bandindices[i][1]+1]:
      temp = l.split()
      if (int(temp[1]) == atomnumber):
	weights.dtot.append(float(temp[orbitalmatching['dtot']]))
	weights.dz2.append(float(temp[orbitalmatching['dz2']]))
	weights.dx2y2.append(float(temp[orbitalmatching['dx2y2']]))
	weights.dxy.append(float(temp[orbitalmatching['dxy']]))
	weights.dxz.append(float(temp[orbitalmatching['dxz']]))
	weights.dyz.append(float(temp[orbitalmatching['dyz']]))
    weightdata.append(weights)
  return weightdata
  
class band:
  def __init__(self, bn):
    self.n = bn
    self.x = []
    self.y = []
    self.z = []
    self.e = []
  
def getbanddata(filename, interestingbands=[1]):
  filehandle = open(filename, 'r')
  lines = filehandle.readlines()
  filehandle.close()

  bandstarts = []
  for i in range(len(lines)):
    if('bandindex:' in lines[i]):
      bandstarts.append(i)

  banddata = []
  for i in interestingbands:
    b = band(i)
    for line in lines[bandstarts[i]+1:bandstarts[i+1]]:
      linedata = line.split()
      b.x.append(float(linedata[0]))
      b.y.append(float(linedata[1]))
      b.z.append(float(linedata[2]))
      b.e.append(float(linedata[4]))
    banddata.append(b)
  return banddata
  
def isfermicrossed(val1, val2):
  iscrossed = 'false'
  if(((val1*val2)<=0)): #if product is smaller than zero sign must have changed
    iscrossed = 'true'
  return iscrossed

def getfermilocation(vec1, vec2, val1, val2):
  fl = []
  if(vec1[0] == vec2[0]):
    x = vec1[0]
    m = (val1 - val2)/(vec1[1] - vec2[1])
    y = vec1[1] - val1/m
    fl = [x,y]
  elif(vec1[1] == vec2[1]):
    m = (val1 - val2)/(vec1[0] - vec2[0])
    x = vec1[0] - val1/m 
    y = vec1[1]
    fl = [x,y]
  else:
    print 'Error. Compairing wrong locations.'
  return fl
  
def getinterpolatedweights(vec1, vec2, vecf, we1, we2):
  w = 0
  if(vec1[0] == vec2[0]):
    m = (we1 - we2)/(vec1[1] - vec2[1])
    w = m*(vecf[1] - vec1[1]) + we1
  elif(vec1[1] == vec2[1]):
    m = (we1 - we2)/(vec1[0] - vec2[0])
    w = m*(vecf[0] - vec1[0]) + we1
  else:
    print 'Error. Compairing wrong locations.'
  return w
  
def getfs(bands, weights):
  n = int(math.sqrt(len(bands[0].x)))
  
  fsbands = []
  for band,weight in zip(bands,weights):
    x = band.x
    y = band.y
    v = band.e
    fs = []
    for i in range(n):
      for j in range(n-1):
	k = j+1
        #search in x-direction
        idx1 = i*n+j
        idx2 = i*n+k
        iscrossed = isfermicrossed(v[idx1], v[idx2])
        if('true' == iscrossed):
	  fl = getfermilocation([x[idx1], y[idx1]], [x[idx2], y[idx2]], v[idx1], v[idx2])
	  we1 = np.array([weight.dz2[idx1],weight.dx2y2[idx1],weight.dxy[idx1],weight.dxz[idx1],weight.dyz[idx1]])
	  we2 = np.array([weight.dz2[idx2],weight.dx2y2[idx2],weight.dxy[idx2],weight.dxz[idx2],weight.dyz[idx2]])
	  w = getinterpolatedweights([x[idx1], y[idx1]], [x[idx2], y[idx2]], fl, we1, we2)
	  maxweightidx = np.argmax(w)
	  maxweight = w[maxweightidx] 
	  fl.extend((maxweightidx+1,maxweight)) #+1 is for convenient plotting in pyxplot using palette option
	  fs.append(fl)
	#search in y-direction
	idx1 = j*n+i
	idx2 = k*n+i
	iscrossed = isfermicrossed(v[idx1], v[idx2])
	if('true' == iscrossed):
	  fl = getfermilocation([x[idx1], y[idx1]], [x[idx2], y[idx2]], v[idx1], v[idx2])
	  we1 = np.array([weight.dz2[idx1],weight.dx2y2[idx1],weight.dxy[idx1],weight.dxz[idx1],weight.dyz[idx1]])
	  we2 = np.array([weight.dz2[idx2],weight.dx2y2[idx2],weight.dxy[idx2],weight.dxz[idx2],weight.dyz[idx2]])
	  w = getinterpolatedweights([x[idx1], y[idx1]], [x[idx2], y[idx2]], fl, we1, we2)
	  maxweightidx = np.argmax(w)
	  maxweight = w[maxweightidx]
	  fl.extend((maxweightidx+1,maxweight)) #+1 is for convenient plotting in pyxplot using palette option
	  fs.append(fl)
    fsbands.append(fs)
  return fsbands
  
def fstofile(filename, fsbands):
  outfilename = filename
  filehandle = open(outfilename, 'w')
  
  filehandle.write('#Fermi Surface points\n#x,y,maxweightidx,maxweight\n#weights are 0=dz2, 1=dx2y2, 2=dxy, 3=dxz, 4=dyz\n')
  for fs in fsbands:
    for vec in fs:
      filehandle.write('%f %f %i %f\n' % (vec[0], vec[1], vec[2], vec[3]))
    filehandle.write('\n\n')
  filehandle.close()
  
def main():
  an = 2 #number of the interesting atom in qtl file
  interestingbands = [30,31,32,33,34,35]
  qtlfilename = 'ggazp.qtl'
  spaghettifilename = 'ggazp.spaghetti_ene'
  outfilename = 'fs.dat'
  orbitalmatching = {'tot': 2, 'stot': 3, 'ptot': 4, 'px': 5, 'py': 6, 'pz': 7, 'dtot': 8, 'dz2': 9, 'dx2y2': 10, 'dxy': 11, 'dxz': 12, 'dyz': 13, 'ftot': 14} #orbital-index matching for ISPLIT=8 
  
  bands= getbanddata(spaghettifilename, interestingbands=interestingbands)
  weights = getqtldata(qtlfilename, orbitalmatching, atomnumber=an, interestingbands=interestingbands)
  fs = getfs(bands, weights)
  fstofile(outfilename, fs)
  
main()
