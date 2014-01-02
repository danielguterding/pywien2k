import os
import numpy as np

def rodriguez(x, y, t):
  #Rodriguez formula, rotate by t around unit vector x starting from y
  x = np.array(x)
  x /= np.linalg.norm(x)
  
  y = np.array(y)
  y /= np.linalg.norm(y)
  
  m = np.zeros((3,3))
  m[0][0] = np.cos(t) + x[0]**2*(1-np.cos(t))
  m[0][1] = x[0]*x[1]*(1-np.cos(t)) - x[2]*np.sin(t)
  m[0][2] = x[1]*np.sin(t) + x[0]*x[2]*(1-np.cos(t))
  m[1][0] = x[2]*np.sin(t) + x[0]*x[1]*(1-np.cos(t))
  m[1][1] = np.cos(t) + x[1]**2*(1-np.cos(t))
  m[1][2] = -x[0]*np.sin(t) + x[1]*x[2]*(1-np.cos(t))
  m[2][0] = -x[1]*np.sin(t) + x[0]*x[2]*(1-np.cos(t))
  m[2][1] = x[0]*np.sin(t) + x[1]*x[2]*(1-np.cos(t))
  m[2][2] = np.cos(t) + x[2]**2*(1-np.cos(t))
  
  return m.dot(y)
  
def main():
  #This script calculates real space directions for the quantization axis of dft calculation with
  #spin orbit coupling. The purpose is to scan the total energy as a function of the quantization axis
  #direction to find the easy axis.
  nang = 37
  angles = np.linspace(0, np.pi, num=nang, endpoint=True)
  
  abh = [rodriguez([0, 0,1], [1,0,0], t) for t in angles] #ab honeycomb plane
  ac  = [rodriguez([0,-1,0], [1,0,0], t) for t in angles] #ac plane
  bc  = [rodriguez([1, 0,0], [0,1,0], t) for t in angles] #bc plane
  
  outfilenames = ['angles_abh.dat', 'angles_ac.dat', 'angles_bc.dat']
  for filename, directions in zip(outfilenames, [abh, ac, bc]):
    filehandle = open(filename, 'w')
    filehandle.write('#deg xdir ydir zdir\n')
    for ang, vec in zip(angles, directions):
      filehandle.write('%f %f %f %f\n' % (ang/np.pi*180, vec[0], vec[1], vec[2]))
    filehandle.close()
  
main()
