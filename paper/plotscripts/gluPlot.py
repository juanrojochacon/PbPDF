#! /usr/bin/python
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

from operator import add, sub
import sys, os, math
import numpy as np

import lhapdf

######################################### Settings ###################################

# Color information
icol = 1
colours = ['r', 'b', 'g']

# Fetch path for replica files  
lhapdf.pathsPrepend("../../fits/NNPbF_REPS09_NNPDF_5_3_replica_2k/")

# Choice of PDF
ipdf = 21
pdfname = "gluon"
outfile = "gluonplot.pdf"

# PDF sets
pset = lhapdf.getPDFSet("NNPbF_REPS09_NNPDF_5_3_replica_2k")
bset = lhapdf.getPDFSet("NNPDF30_nlo_as_0118")


####################################FIGURE SETUP##################################

 # Setup gridspec
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 0.5])
gs.update(wspace=0.00, hspace=0.00)

# Setup figure
fig = plt.figure()
lax =  fig.add_subplot(gs[0])
ax = fig.add_subplot(gs[1])

rlax =  fig.add_subplot(gs[2])
rax = fig.add_subplot(gs[3])


# Gridlines
ax.xaxis.grid(True)
ax.yaxis.grid(True)

# Gridlines
lax.xaxis.grid(True)
lax.yaxis.grid(True)

# Gridlines
rax.xaxis.grid(True)
rax.yaxis.grid(True)

# Gridlines
rlax.xaxis.grid(True)
rlax.yaxis.grid(True)

# Axis formatting
plt.setp(ax.get_yticklabels(), visible=False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(lax.get_xticklabels(), visible=False)

lax.set_xscale('log')

ax.xaxis.set_major_locator(MaxNLocator(5,prune='lower'))
lax.yaxis.set_major_locator(MaxNLocator(8,prune='lower'))
ax.yaxis.set_major_locator(MaxNLocator(8,prune='lower'))

lax.set_ylabel("f_" + pdfname)
ax.set_xlabel("x")
lax.set_xlabel("x")

# Axis formatting
plt.setp(rax.get_yticklabels(), visible=False)
rlax.set_xscale('log')

rax.xaxis.set_major_locator(MaxNLocator(5,prune='lower'))

rlax.set_ylabel("R_" + pdfname)
rax.set_xlabel("x")
rlax.set_xlabel("x")

########################################### Plot data ###########################################

# Kinematics
xs = [x for x in np.logspace(-7, -0.0001, 500)]
q = 1.3

yvals = []
rvals = []
bvals = []

nreps=pset.size
interval = (nreps-0.68*nreps)/2

for i in xrange(0,pset.size):
  prep = pset.mkPDF(i)
  brep = bset.mkPDF(i+1)

  repyvals = [ prep.xfxQ(ipdf, x, q) for x in xs ]
  repbvals = [ brep.xfxQ(ipdf, x, q) for x in xs ]
  reprvals = [ prep.xfxQ(ipdf, x, q)/brep.xfxQ(ipdf, x, q) for x in xs ]

  yvals.append(repyvals)
  rvals.append(reprvals)
  bvals.append(repbvals)

  ax.plot(xs,repyvals, color = colours[icol], alpha=0.05)
  lax.plot(xs,repyvals, color = colours[icol], alpha=0.05)
  rax.plot(xs,reprvals, color = colours[icol], alpha=0.05)
  rlax.plot(xs,reprvals, color = colours[icol], alpha=0.05)

# Compute averages
yCV = []; yER = []; 
y84 = []; y16 = []; 
rCV = []; rER = [];
r84 = []; r16 = [];
bCV = []; bER = [];
rbCV = []; rbER = [];
yvals_trans = map(list, zip(*yvals))
rvals_trans = map(list, zip(*rvals))
bvals_trans = map(list, zip(*bvals))

for replicas in yvals_trans:
  srtreps = np.sort(replicas)
  yCV.append(np.mean(srtreps))
  yER.append(np.std(srtreps))
  y84.append(srtreps[nreps-interval])
  y16.append(srtreps[interval])

for replicas in rvals_trans:
  srtreps = np.sort(replicas)
  rCV.append(np.mean(srtreps))
  rER.append(np.std(srtreps))
  r84.append(srtreps[nreps-interval])
  r16.append(srtreps[interval])

rbCV = [ 1 for x in xs ]
for replicas in bvals_trans:
  srtreps = np.sort(replicas)
  bCV.append(np.mean(srtreps))
  bER.append(np.std(srtreps))

CVup = map(add, yCV, yER); rCVup = map(add, rCV, rER)
CVdn = map(sub, yCV, yER); rCVdn = map(sub, rCV, rER)
bCVup = map(add, bCV, bER); bCVdn = map(sub, bCV, bER)

# PROBABLY NOT CORRECT - NEED TO DIVIDE BY CV
rbER = [a/b for a,b in zip(bER,bCV)]
rbCVup = map(add, rbCV, rbER); rbCVdn = map(sub, rbCV, rbER)

########################################### Plotting ###########################################


# Central values
ax.plot(xs,yCV, color = colours[icol], alpha=0.8, label = "NNPbF")
lax.plot(xs,yCV, color = colours[icol], alpha=0.8, label = "NNPbF")
rax.plot(xs,rCV, color = colours[icol], alpha=0.8, label = "NNPbF")
rlax.plot(xs,rCV, color = colours[icol], alpha=0.8, label = "NNPbF")

# error bars
ax.fill_between(xs, CVdn, CVup, alpha=0.2,
                      facecolor=colours[icol], linewidth=0, color=colours[icol])

rax.fill_between(xs, rCVdn, rCVup, alpha=0.2,
                      facecolor=colours[icol], linewidth=0, color=colours[icol])

ax.fill_between(xs, bCVdn, bCVup, alpha=0.2,
                      facecolor=colours[icol+1], linewidth=0, color=colours[icol])
rax.fill_between(xs, rbCVdn, rbCVup, alpha=0.2,
                      facecolor=colours[icol+1], linewidth=0, color=colours[icol])



# error bars
lax.fill_between(xs, CVdn, CVup, alpha=0.2,
                      facecolor=colours[icol], linewidth=0, color=colours[icol])

lax.fill_between(xs, bCVdn, bCVup, alpha=0.2,
                      facecolor=colours[icol+1], linewidth=0, color=colours[icol])

rlax.fill_between(xs, rCVdn, rCVup, alpha=0.2,
                      facecolor=colours[icol], linewidth=0, color=colours[icol])

rlax.fill_between(xs, rbCVdn, rbCVup, alpha=0.2,
                      facecolor=colours[icol+1], linewidth=0, color=colours[icol])

# 68CL
ax.plot(xs,y84, color = colours[icol], alpha=0.8, linestyle='-.')
lax.plot(xs,y84, color = colours[icol], alpha=0.8, linestyle='-.')
rax.plot(xs,r84, color = colours[icol], alpha=0.8, linestyle='-.')
rlax.plot(xs,r84, color = colours[icol], alpha=0.8, linestyle='-.')

  # 68CL
ax.plot(xs,y16, color = colours[icol], alpha=0.8, linestyle='-.')
lax.plot(xs,y16, color = colours[icol], alpha=0.8, label = "68% C.I", linestyle='-.')

rax.plot(xs,r16, color = colours[icol], alpha=0.8, linestyle='-.')
rlax.plot(xs,r16, color = colours[icol], alpha=0.8, label = "68% C.I", linestyle='-.')



########################################### Finalize ###########################################

# set limits
rax.set_xlim([0.1, 1])
rlax.set_xlim([1E-5,0.1])

# set limits
ax.set_xlim([0.1, 1])
lax.set_xlim([1E-5,0.1])

rax.set_ylim( [0,2])
rlax.set_ylim([0,2])

ax.set_ylim( [-1,3])
lax.set_ylim([-1,3])

fig.savefig(outfile)




