import sys
import os

oldname = sys.argv[1]
newname = sys.argv[2]
for f in os.listdir('.'):
  if f.startswith(oldname):
    os.rename(f, f.replace(oldname, newname))
