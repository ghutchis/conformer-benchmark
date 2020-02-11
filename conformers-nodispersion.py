 #!/usr/bin/env python

from __future__ import print_function

import sys, os
import glob
import time

# use open babel 3.0
from openbabel import pybel

# data with per-atom SCF energies
import atomization

hartree_to_kcal = 627.509469
eV_to_kcal = 23.06035

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def key(dir, file):
    subdir = dir
    if '/' in dir:
        # remove the "CHG" or "Neutral" top-level directory
        # if it's part of the string
        top, subdir = dir.split('/')
    if '.' in file:
        # remove the extension if it's still part of the filename
        file = file.split('.')[0]
    # returns something like astex_1gm8-rmsd112-opt
    return "{}-{}".format(subdir, file)

# we need to read through all the energy files

ccsdt = {}
with open('energies/ccsdt.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) > 7:
            continue
        k = key(items[0], items[1])
        ccsdt[k] = float(items[6]) # in Hartrees

pbe = {}
with open('energies/pbe-nod.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 Total Energy       :        -1502.
        items = line.split()
        if len(items) < 7:
            continue
        k = key(items[0], items[1])
        pbe[k] = float(items[5]) # in Hartrees

pbe_svp = {}
with open('energies/pbe-svp-nod.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 Total Energy       :        -1502.
        items = line.split()
        if len(items) < 7:
            continue
        k = key(items[0], items[1])
        pbe_svp[k] = float(items[5]) # in Hartrees

wb97 = {}
with open('energies/wb97xd-nod.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) != 9:
            continue
        k = key(items[0], items[1])
        wb97[k] = float(items[5]) # in Hartrees

b3lyp_tz = {}
with open('energies/b3lyp-tz-nod.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) < 7:
            continue
        k = key(items[0], items[1])
        b3lyp_tz[k] = float(items[5]) # in Hartrees

b3lyp_svp = {}
with open('energies/b3lyp-svp-nod.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) < 7:
            continue
        k = key(items[0], items[1])
        b3lyp_svp[k] = float(items[5]) # in Hartrees

# okay, now we use pybel to open and parse the .mol files
print("name,geom,natoms,dlpno,wb97,pbe,pbeSVP,b3lypTZ,b3lypSVP")
for filename in glob.iglob("/Users/ghutchis/conf/*_jobs/*/*opt.mol"):
    name = filename.split('/')[-2]
    geom = filename.split('/')[-1]

    mol = next(pybel.readfile("mol", filename))
    elements = [0] * 118
    for atom in mol.atoms: # how many of each element are there?
        elements[atom.atomicnum] = elements[atom.atomicnum] + 1

    # okay, we can use elements[] to calculate atomization energies if needed
    dlpnoE = mp2E = b3lypSVPE = b3lypTZE = float('nan')
    wb97E = pbeE = gfnE = gfn2E = pm7E = pm7HOFE = float('nan')
    pbeSVPE = pbeh3cE = b973cE = float('nan')
    ani1E = ani1ccE = ani2E = float('nan')
    k = key(name, geom)

    # now each of the energies / dictionaries in turn
    if k in ccsdt:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.dlpnoE[e]
        dlpnoE = (totalAtomE - float(ccsdt[k])) * hartree_to_kcal
    if k in pbe:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.pbeE[e]
        pbeE = (totalAtomE - float(pbe[k])) * hartree_to_kcal
    if k in pbe_svp:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.pbe_svpE[e]
        pbeSVPE = (totalAtomE - float(pbe_svp[k])) * hartree_to_kcal
    if k in wb97:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.wb97E[e]
        wb97E = (totalAtomE - float(wb97[k])) * hartree_to_kcal
    if k in b3lyp_tz:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.b3lypTZE[e]
        b3lypTZE = (totalAtomE - float(b3lyp_tz[k])) * hartree_to_kcal
    if k in b3lyp_svp:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.b3lypE[e]
        b3lypSVPE = (totalAtomE - float(b3lyp_svp[k])) * hartree_to_kcal

    # energies
    print(name, geom, len(mol.atoms), dlpnoE, wb97E, pbeE, pbeSVPE, b3lypTZE, b3lypSVPE, sep=',')
