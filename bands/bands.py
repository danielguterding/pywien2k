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
    for l in lines[bandindices[i-1][0]:bandindices[i-1][1]+1]:
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
    self.k = []
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
    partlines = []
    if(i < (len(bandstarts))):
      partlines = lines[bandstarts[i-1]+1:bandstarts[i]]
    elif(i == (len(bandstarts))):
      partlines = lines[bandstarts[i-1]+1:]
    else:
      print 'Error. Band index out of available range.'
      exit()
    for line in partlines:
      linedata = line.split()
      b.x.append(float(linedata[0]))
      b.y.append(float(linedata[1]))
      b.z.append(float(linedata[2]))
      b.k.append(float(linedata[3]))
      b.e.append(float(linedata[4]))
    banddata.append(b)
  return banddata
  
def writebanddata(outfilename, banddata, weights):
  of = open(outfilename, 'w')
  of.write('#band index, kindex, energy, weight dtot, dz2, dx2y2, dxy, dxz, dyz\n')
  for band,weight in zip(banddata, weights):
    for i in range(len(band.k)):
      of.write('%i %f %f %f %f %f %f %f %f\n' % (band.n, band.k[i], band.e[i], weight.dtot[i], weight.dz2[i], weight.dx2y2[i], weight.dxy[i], weight.dxz[i], weight.dyz[i]))
    of.write('\n\n')
  of.close()

def main():
  an = 1 #number of the interesting atom in qtl file
  interestingbands = [i for i in range(5,18)]
  
  qtlfilename = 'NiO.qtl'
  spaghettifilename = 'NiO.spaghetti_ene'
  bandsfilename = 'NiO.bands'
  orbitalmatching = {'tot': 2, 'stot': 3, 'ptot': 4, 'px': 5, 'py': 6, 'pz': 7, 'dtot': 8, 'dz2': 9, 'dx2y2': 10, 'dxy': 11, 'dxz': 12, 'dyz': 13, 'ftot': 14} #orbital-index matching for ISPLIT=8 
  
  bands = getbanddata(spaghettifilename, interestingbands=interestingbands)
  weights = getqtldata(qtlfilename, orbitalmatching, atomnumber=an, interestingbands=interestingbands)
  writebanddata(bandsfilename, bands, weights)
  
main()