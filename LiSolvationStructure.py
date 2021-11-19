#!/usr/bin/python
#--------------------------------------------------------------------------------------------
# code used to analyze solvation environment of lithium ions in different electrolytes
# for the paper "Suspension electrolyte with modified Li+ solvation environment for lithium metal batteries"
# the default trajectory and MD setting are run.xtc and run.tpr, respectively
# For questions, reach Paul Rudnicki at prudnick@stanford.edu
#--------------------------------------------------------------------------------------------


import math as m
import os
import numpy as np
import argparse

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TPR, XTC

#Natural Constants
kB=1.38064852e-23
e=1.60217662e-19
T=300.0
beta=1.0/(kB*T)
dr=0.002

#process inputs
parser = argparse.ArgumentParser()
parser.add_argument("-ani",help="anion", default='TFS')
parser.add_argument("-sol",help="solvent", default='EC')
parser.add_argument("-b",help="begin time (ps)", default=0,type=float)
parser.add_argument("-dt",help="time step in trajectory (ps)", default=1,type=float)
parser.add_argument("-dta",help="time step for analysis (ps)", default=10,type=float)
args = parser.parse_args()

ani=args.ani
sol=args.sol
beg=args.b
dt=args.dt
dta=args.dta

sols=[ani,sol]

if sol=='ED':
   sols=[ani,'EC','DEC']
if sol=='EDF':
   sols=[ani,'EC','DEC','FEC']

atmani="O"
if ani=='PF6':
    atmani="P"
    #atmani="F"
if ani=='BF4':
    atmani="F"
if ani=='I':
    atmani="I"
if ani=='FSI':
    atmani="O"

nsols=len(sols)

maxSol=100
maxAni=100

if nsols==4:
   counts=np.zeros([maxAni,maxSol,maxSol,maxSol])
elif nsols==3:
   counts=np.zeros([maxAni,maxSol,maxSol])
elif nsols==2:
   counts=np.zeros([maxAni,maxSol])
else:
   print("unknown number of solvents, exiting...")
   exit()

#--------------------------------------------------------------
# setting cutoff value, in A, manually, for testing
#--------------------------------------------------------------

#rcuts=3.0*np.ones(nsols)

#--------------------------------------------------------------
# determing cutoff value from RDFs, used in the manuscript
#--------------------------------------------------------------

rcuts=np.zeros(nsols)
frdfLiN='rdf.LI.{:s}of{:s}.xvg'.format(atmani,ani)
dat0=np.fromfile(frdfLiN,sep=" ")
print(dat0)
print( len(dat0) )
rdfLiN=dat0.reshape(int(len(dat0)/2),2)
iMax=np.argmax(rdfLiN[0:int(0.6/dr),1])
iMin=np.argmin(rdfLiN[iMax:int(0.6/dr),1])
rcuts[0]=10*rdfLiN[iMax+iMin,0]

for isol in range(1,nsols):
   si=sols[isol]
   frdfLiOE='rdf.LI.OEof'+si+'.xvg'
   frdfLiOC='rdf.LI.OCof'+si+'.xvg'
   if os.path.isfile(frdfLiOC):
      frdfLiO=frdfLiOC
   else:
      frdfLiO=frdfLiOE

   dat0=np.fromfile(frdfLiO,sep=" ")
   rdfLiO=dat0.reshape(int(len(dat0)/2),2)
   iMax=np.argmax(rdfLiO[:,1])
   iMin=np.argmin(rdfLiO[iMax:int(0.6/dr),1])
   rcuts[isol]=10*rdfLiO[iMax+iMin,0]

print(rcuts)

u = mda.Universe("run.tpr", "run.xtc", in_memory=False)
li = u.select_atoms("resname LI")

nFrame = len(u.trajectory)

#print(nFrame)
#print (resIDSol)
print(li.atoms)

matRes=np.zeros((nFrame,1))

for ts in u.trajectory[int(beg/dt):]:
   if (u.trajectory.time > beg) and (u.trajectory.time %dta==0 ):
#      print("Frame: {0:5d}, Time: {1:8.3f} ps, nsol: {2:5d} ".format(ts.frame, u.trajectory.time, 0))
      for a in li.atoms:
         agLi1=mda.core.groups.AtomGroup([a])
         nSolRes=np.zeros((nsols,), dtype=int)
         nSolRes2=np.zeros((nsols,), dtype=int)
         for isol in range(0,nsols):
            si=sols[isol]
            solAtms = u.select_atoms("resname {0:5s} and (not name C* B* H* S* F*) and around {1:8.3f} group li1".format(si,rcuts[isol]),li1=agLi1)
            solRes = solAtms.residues
            resIDSol = solRes.resids
            if isol==0:
                #print(resIDSol)
                print(solRes)
            solAtms2= u.select_atoms("resname {0:5s} and name O and around {1:8.3f} group li1".format(si,rcuts[isol]),li1=agLi1)
            solRes2= solAtms2.residues
            resIDSol2 = solRes2.resids
            #print(solAtms)
            #print(solRes)
            #print(resIDSol)
            nSolRes[isol] = len(resIDSol)
            nSolRes2[isol] = len(resIDSol2)
            #print("Frame: {0:5d}, Time: {1:8.3f} ps, nsol: {2:5d} ".format(ts.frame, u.trajectory.time, nSolRes[isol]))
         if nsols==4:
            counts[nSolRes[0],nSolRes[1],nSolRes[2],nSolRes[3]] += 1.0
         elif nsols==3:
            counts[nSolRes[0],nSolRes[1],nSolRes[2]] += 1.0
         elif nsols==2:
            counts[nSolRes[0],nSolRes[1]] += 1.0

sumCounts=np.sum(counts)
freqs=counts/sumCounts*100

fo=open('coordNumFreq.log',"w")

if nsols==4:
 for i1 in range(maxAni):
  for i2 in range(maxSol):
   for i3 in range(maxSol):
    for i4 in range(maxSol):
     if freqs[i1,i2,i3,i4]>=0.1:
      print("{:5d} {:5d} {:5d} {:5d} {:8.3f}".format(i1,i2,i3,i4,freqs[i1,i2,i3,i4]) )
      fo.write("{:5d} {:5d} {:5d} {:5d} {:8.3f} \n".format(i1,i2,i3,i4,freqs[i1,i2,i3,i4]) )
elif nsols==3:
 for i1 in range(maxAni):
  for i2 in range(maxSol):
   for i3 in range(maxSol):
     if freqs[i1,i2,i3]>=0.1:
      print("{:5d} {:5d} {:5d} {:8.3f}".format(i1,i2,i3,freqs[i1,i2,i3]) )
      fo.write("{:5d} {:5d} {:5d} {:8.3f} \n".format(i1,i2,i3,freqs[i1,i2,i3]))
elif nsols==2:
 for i1 in range(maxAni):
  for i2 in range(maxSol):
     if freqs[i1,i2]>=0.1:
      print("{:5d} {:5d} {:8.3f}".format(i1,i2,freqs[i1,i2]) )
      fo.write("{:5d} {:5d} {:8.3f} \n".format(i1,i2,freqs[i1,i2]))
fo.close()
