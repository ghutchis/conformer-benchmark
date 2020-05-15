 #!/usr/bin/env python

from __future__ import print_function

import sys
import os
import glob
import time

from statistics import mean, median
from sklearn.utils import resample
import numpy as np

# bootstrap confidence intervals
# repeatedly resample the data into random subsets
# .. take the mean and medians of the subset
# .. you can then use numpy.percentile() to get CI
def boot_sample(list, iterations=1000, fraction=0.75):
    size = int( len(list) * fraction)
    means = []
    medians = []
    for i in range(iterations):
        subset = resample(list, n_samples=size)
        means.append(mean(subset))
        medians.append(median(subset))
    return means, medians

# uses open babel 3.0
from openbabel import pybel
mmff = pybel._forcefields["mmff94"]
mmffS = pybel._forcefields["mmff94s"] # unused
uff = pybel._forcefields["uff"]
gaff = pybel._forcefields["gaff"]

# open up some timing files:
mmffTime = open("timing/mmff-batch.csv", 'w')
mmff_times = []
uffTime = open("timing/uff-batch.csv", 'w')
uff_times = []
gaffTime = open("timing/gaff-batch.csv", 'w')
gaff_times = []

kj_to_kcal = 1.0 / 4.184

# evaluate a batch of geometries with force fields
current_mol_name = ""
for filename in glob.iglob("geometries/*_jobs/*/*opt.mol"):
    name = filename.split('/')[-2]
    geom = filename.split('/')[-1]

    # we read in a new file
    mol = next(pybel.readfile("mol", filename))
    if name != current_mol_name:
        # we have to set up the force field parameters
        mmff.Setup(mol.OBMol)
        uff.Setup(mol.OBMol)
        gaff.Setup(mol.OBMol)
        current_mol_name = name

    # if it's the same molecule, we avoid Setup()
    # but need to update coordinates (which should be timed)

    # get the force field energies
    t0 = time.perf_counter()
    mmff.SetCoordinates(mol.OBMol)
    # false = don't calculate gradients
    mmffE = mmff.Energy(False) # in kcal/mol
    t1 = time.perf_counter()
    print(name, geom, t1-t0, sep=',', file=mmffTime)
    mmff_times.append(t1 - t0)

    t0 = time.perf_counter()
    uff.SetCoordinates(mol.OBMol)
    # false = don't calculate gradients
    uffE = uff.Energy(False) * kj_to_kcal
    t1 = time.perf_counter()
    print(name, geom, t1-t0, sep=',', file=uffTime)
    uff_times.append(t1-t0)

    t0 = time.perf_counter()
    gaff.SetCoordinates(mol.OBMol)
    # false = don't calculate gradients
    gaffE = gaff.Energy(False) * kj_to_kcal
    t1 = time.perf_counter()
    print(name, geom, t1-t0, sep=',', file=gaffTime)
    gaff_times.append(t1-t0)

mmffTime.close()
uffTime.close()
gaffTime.close()

# okay, now lets report medians and CI for each
means, medians = boot_sample(mmff_times)
med = median(mmff_times)
low_cim = med - np.percentile(medians, 2.5)
high_cim = np.percentile(medians, 97.5) - med
print("MMFF94: ", med, low_cim, high_cim)

means, medians = boot_sample(gaff_times)
med = median(gaff_times)
low_cim = med - np.percentile(medians, 2.5)
high_cim = np.percentile(medians, 97.5) - med
print("GAFF: ", med, low_cim, high_cim)

means, medians = boot_sample(uff_times)
med = median(uff_times)
low_cim = med - np.percentile(medians, 2.5)
high_cim = np.percentile(medians, 97.5) - med
print("UFF: ", med, low_cim, high_cim)
