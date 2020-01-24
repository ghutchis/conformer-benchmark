 #!/usr/bin/env python

from __future__ import print_function

import sys, os
import glob
import time

import pybel
mmff = pybel._forcefields["mmff94"]
mmffS = pybel._forcefields["mmff94s"]
uff = pybel._forcefields["uff"]
gaff = pybel._forcefields["gaff"]

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

mp2 = {}
with open('energies/mp2.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) > 7:
            continue
        k = key(items[0], items[1])
        mp2[k] = float(items[6]) # in Hartrees

pbe = {}
with open('energies/pbe.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) > 7:
            continue
        k = key(items[0], items[1])
        pbe[k] = float(items[6]) # in Hartrees

pbe_svp = {}
with open('energies/pbe-svp.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) > 7:
            continue
        k = key(items[0], items[1])
        pbe_svp[k] = float(items[6]) # in Hartrees

wb97 = {}
with open('energies/wb97xd.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) != 7:
            continue
        k = key(items[0], items[1])
        wb97[k] = float(items[6]) # in Hartrees

b3lyp_tz = {}
with open('energies/b3lyp-tz.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) > 7:
            continue
        k = key(items[0], items[1])
        b3lyp_tz[k] = float(items[6]) # in Hartrees

b3lyp_svp = {}
with open('energies/b3lyp-svp.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) > 7:
            continue
        k = key(items[0], items[1])
        b3lyp_svp[k] = float(items[6]) # in Hartrees

b973c = {}
with open('energies/b97-3c.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) > 7:
            continue
        k = key(items[0], items[1])
        b973c[k] = float(items[6]) # in Hartrees

pbeh3c = {}
with open('energies/pbeh-3c.txt') as file:
    for line in file:
        # astex_1gm8 rmsd112-opt.out.bz2 FINAL SINGLE POINT ENERGY     -1501.5
        items = line.split()
        if len(items) > 7:
            continue
        k = key(items[0], items[1])
        pbeh3c[k] = float(items[6]) # in Hartrees

pm7 = {} # total energies to convert to atomization energy
with open('energies/pm7-total.txt') as file:
    for line in file:
        items = line.split()
        if len(items) != 6:
            continue

        # CHG_jobs/astex_1gm8/rmsd1-mmff.pm7:
        dir = items[0].split('/')[1]
        file = items[0].split('/')[2]
        k = key(dir,file)
        pm7[k] = float(items[4]) # in eV

pm7hof = {} # heats of formation
with open('energies/pm7-hof.txt') as file:
    for line in file:
        items = line.split()
        if len(items) != 11:
            continue

        # CHG_jobs/astex_1gm8/rmsd1-mmff.pm7:
        dir = items[0].split('/')[1]
        file = items[0].split('/')[2]
        k = key(dir,file)
        pm7hof[k] = float(items[6]) # in kcal/mol

gfn0 = {}
with open('energies/gfn0.txt') as file:
    for line in file:
        # parse it
        items = line.split()
        k = key(items[0], items[1])
        gfn0[k] = float(items[5]) # in Hartrees

# older XTB fortunately outputs atomization energies in Hartree and kcal/mol
gfn = {}
with open('energies/gfn.txt') as file:
    for line in file:
        # parse it
        items = line.split()
        k = key(items[0], items[1])
        gfn[k] = float(items[5]) # in kcal/mol

gfn2 = {}
with open('energies/gfn2.txt') as file:
    for line in file:
        # parse it
        items = line.split()
        k = key(items[0], items[1])
        gfn2[k] = float(items[5]) # in kcal/mol

ani1 = {}
with open('energies/ani-1x.txt') as file:
    for line in file:
        items = line.split(',')
        name = items[0].split('/')[-1]
        dir = '/'.join(items[0].split('/')[0:2])
        k = key(dir, name)
        ani1[k] = float(items[1]) # in kcal/mol

ani1cc = {}
with open('energies/ani-1ccx.txt') as file:
    for line in file:
        items = line.split(',')
        name = items[0].split('/')[-1]
        dir = '/'.join(items[0].split('/')[0:2])
        k = key(dir, name)
        ani1cc[k] = float(items[1]) # in kcal/mol

ani2 = {}
with open('energies/ani-2x.txt') as file:
    for line in file:
        items = line.split(',')
        name = items[0].split('/')[-1]
        dir = '/'.join(items[0].split('/')[0:2])
        k = key(dir, name)
        ani2[k] = float(items[1]) # in kcal/mol


# open up some timing files:
mmffTime = open("timing/mmff-time.csv", 'w')
uffTime = open("timing/uff-time.csv", 'w')
gaffTime = open("timing/gaff-time.csv", 'w')

# okay, now we use pybel to open and parse the .mol files
print("name,geom,natoms,dlpno,mp2,wb97,b973c,pbe,pbeSVP,pbeh3c,b3lypTZ,b3lypSVP,gfn0,gfn1,gfn2,pm7E,pm7HOF,mmff,mmffNew,uff,gaff,ani1x,ani1cc,ani2")
for filename in glob.iglob("/Users/ghutchis/conf/*_jobs/*/*.mol"):
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
    if k in mp2:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.mp2E[e]
        mp2E = (totalAtomE - float(mp2[k])) * hartree_to_kcal
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
    if k in pbeh3c:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.pbeh3cE[e]
        pbeh3cE = (totalAtomE - float(pbeh3c[k])) * hartree_to_kcal
    if k in wb97:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.wb97E[e]
        wb97E = (totalAtomE - float(wb97[k])) * hartree_to_kcal
    if k in b973c:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.b973cE[e]
        b973cE = (totalAtomE - float(b973c[k])) * hartree_to_kcal
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
    if k in gfn0:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.gfn0E[e]
        gfn0E = (totalAtomE - float(gfn0[k])) * hartree_to_kcal

    if k in pm7:
        totalAtomE = 0.0
        for e in range(len(elements)):
            totalAtomE = totalAtomE + elements[e] * atomization.pm7E[e]
        pm7E = (totalAtomE - float(pm7[k])) * eV_to_kcal

    # already in kcal/mol
    if k in gfn:
        gfnE = gfn[k]
    if k in gfn2:
        gfn2E = gfn2[k]
    if k in ani1:
        ani1E = ani1[k]
    if k in ani1cc:
        ani1ccE = ani1cc[k]
    if k in ani2:
        ani2E = ani2[k]
    if k in pm7hof:
        pm7HOFE = pm7hof[k]

    # get the force field energies
    t0 = time.perf_counter()
    mmff.Setup(mol.OBMol)
    mmffE = mmff.Energy() # in kcal/mol
    t1 = time.perf_counter()
    print(name, geom, t1-t0, sep=',', file=mmffTime)

    t0 = time.perf_counter()
    uff.Setup(mol.OBMol)
    uffE = uff.Energy()
    t1 = time.perf_counter()
    print(name, geom, t1-t0, sep=',', file=uffTime)

    t0 = time.perf_counter()
    gaff.Setup(mol.OBMol)
    gaffE = gaff.Energy()
    t1 = time.perf_counter()
    print(name, geom, t1-t0, sep=',', file=gaffTime)

    # energies
    print(name, geom, len(mol.atoms), dlpnoE, mp2E, wb97E, b973cE, pbeE, pbeSVPE, pbeh3cE, b3lypTZE, b3lypSVPE, gfn0E, gfnE, gfn2E, pm7E, pm7HOFE, mmffE, uffE, gaffE, ani1E,ani1ccE,ani2E, sep=',')


    #    print(items[0], items[1], ccE, items[2], items[3], pm7E, gfnE, gfn2E, items[4].strip(), sep=',')

mmffTime.close()
uffTime.close()
gaffTime.close()
