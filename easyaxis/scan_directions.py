import os
import sys
import subprocess
import numpy as np

class direc:
  def __init__(self):
    self.deg = 0
    self.rad = 0
    self.vec = np.zeros(3)
  def set_angle_by_rad(self,an):
    self.rad = an
    self.deg = an*180/np.pi
  def set_angle_by_deg(self,an):
    self.deg = an
    self.rad = an/180.0*np.pi
  def set_vector(self,x,y,z):
    self.vec = np.array([x,y,z])
    self.x = self.vec[0]
    self.y = self.vec[1]
    self.z = self.vec[2]

def get_directions(filename):
  filehandle = open(filename, 'r')
  lines = filehandle.readlines()[1:]
  filehandle.close()
  
  directions = []
  for l in lines:
    an, x, y, z = l.split()
    d = direc()
    d.set_angle_by_deg(float(an))
    d.set_vector(float(x), float(y), float(z))
    directions.append(d)
  return directions
  
def generate_insofile(materialname,h,k,l):
  filename = materialname + '.inso'
  filehandle = open(filename, 'w')
  filehandle.write("""WFFIL
4  0  0                 llmax,ipr,kpot
-10  1.5                Emin, Emax
%1.3f %1.3f %1.3f      h,k,l (direction of magnetization)
 0                       number of atoms with RLO
0 0      number of atoms without SO, atomnumbers""" % (h,k,l)
)
  
  filehandle.close()
  
def generate_machines_file(machinename):
  filename = '.machines'
  filehandle = open(filename, 'w')
  filehandle.write("""1:%s
1:%s
1:%s
1:%s
granularity:1
extrafine:1""" % (machinename, machinename, machinename, machinename))
  
def run_wien2k_convergence():
  command = 'runsp_lapw -p -so -orb -ec 0.0001 -cc 0.0001' 
  p = subprocess.Popen(command.split()) 
  p.wait()
  
def get_converged_energy(materialname):
  scffilename = materialname + '.scf'
  command = 'grep :ENE %s | tail -n 1' % (scffilename) #get the last energy value in the scf file
  p = subprocess.Popen(command.split(), shell=False, stdout=subprocess.PIPE) 
  p.wait()
  stdout, stderr = p.communicate()
  energy = float(stdout.split()[-1].strip())
  return energy
  
def get_converged_magnetic_moment(materialname, momentindex):
  scffilename = materialname + '.scf'
  command = 'grep :MMI%03i %s | tail -n 1' % (momentindex, scffilename) #get the value for the magnetic moment in the scf file
  p = subprocess.Popen(command.split(), shell=False, stdout=subprocess.PIPE) 
  p.wait()
  stdout, stderr = p.communicate()
  moment = float(stdout.split()[-1].strip())
  return moment
  
def remove_broyden_files():
  command = 'rm *broyd*'
  p = subprocess.Popen(command.split(), shell=False, stdout=subprocess.PIPE) 
  p.wait()
  
def main():
  angfilename = sys.argv[1]
  materialname = os.getcwd().split('/')[-1] #get the name of the working folder without the path
  directions = get_directions(angfilename)
  
  outfilename = os.path.splitext(angfilename)[0] + '_res.dat'
  outfilehandle = open(outfilename, 'w')
  outfilehandle.write('#deg, xdir, ydir, zdir, energy in meV, magnetic moment 001 in mu bohr\n')
  outfilehandle.close()
  
  firstenergy = None

  for d in directions:
    generate_insofile(materialname, d.x, d.y, d.z)
    generate_machines_file('fullen')
    run_wien2k_convergence()
    energy = get_converged_energy(materialname)
    magneticmoment = get_converged_magnetic_moment(materialname, 1)
    if(None == firstenergy): #save energy of first point so that output can be written relative to that one
      firstenergy = energy
    outfilehandle = open(outfilename, 'a')
    outfilehandle.write('%f %f %f %f %f %f\n' % (d.deg, d.x, d.y, d.z, (energy-firstenergy)*13600, magneticmoment))
    outfilehandle.close()
    remove_broyden_files()
    
main()