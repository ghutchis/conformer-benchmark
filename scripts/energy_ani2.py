# -*- coding: utf-8 -*-
"""
Computing energies using TorchANI
"""

###############################################################################
# To begin with, let's first import the modules we will use:
from __future__ import print_function

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

import torch
import torchani

from timeit import default_timer as timer

import sys
import os

# path to ANI-2x subdirectories
path = '/Users/ghutchis/Devel/torchani/torchani/resources/ani-2x'

# for ANI-2x from sae
atom_selfE = {
'H': -0.5975740895940036,
'C': -38.08963824205658,
'N': -54.71039678680539,
'O': -75.19115551439664,
'S': -398.1787966261814,
'F': -99.79839227891235,
'Cl': -460.206968972609
}

start = timer()

#######
# We're going to read in an XYZ file specified on the command-line
coords = []
self_energy = 0.0
atoms = []
if len(sys.argv) > 1:
    with open(sys.argv[1]) as xyz_file:
        num_atoms = int(xyz_file.readline().rstrip())
        xyz_file.readline() # skip the smiles or title line
        for i in range(num_atoms):
            element, x, y, z = xyz_file.readline().split()
            atoms.append(element)
            coords.append([float(x),float(y),float(z)])
            self_energy += atom_selfE[element]

# optionally check for non-parameterized elements

# Manually specify the device we want TorchANI to run:
device = torch.device('cpu')

###############################################################################
# Let's now load the built-in ANI-1 models. The builtin ANI-1ccx contains 8
# models trained with diffrent initialization. Predicting the energy and force
# using the average of the 8 models outperform using a single model, so it is
# always recommended to use an ensemble, unless the speed of computation is an
# issue in your application.

# (not sure if these need the level shifters...)
# model = torchani.models.ANI1ccx()
# model = torchani.models.ANI1x()


# this whole pile is to load ANI-2x (since it's not public)
const_file = path + '/rHCNOSFCl-5.1R_16-3.5A_a8-4.params'  # noqa: E501
consts = torchani.neurochem.Constants(const_file)
aev_computer = torchani.AEVComputer(**consts)
sae_file = path + '/sae_linfit.dat'  # noqa: E501
energy_shifter = torchani.neurochem.load_sae(sae_file)
model_prefix = path + '/train'  # noqa: E501
ensemble = torchani.neurochem.load_model_ensemble(consts.species, model_prefix, 8)  # noqa: E501

model = torch.nn.Sequential(aev_computer, ensemble, energy_shifter)
# done (that's ANI-2x)

###############################################################################
# Now let's define the coordinate and species. If you just want to compute the
# energy and force for a single structure like in this example, you need to
# make the coordinate tensor has shape ``(1, Na, 3)`` and species has shape
# ``(1, Na)``, where ``Na`` is the number of atoms in the molecule, the
# preceding ``1`` in the shape is here to support batch processing like in
# training. If you have ``N`` different structures to compute, then make it
# ``N``.
coordinates = torch.tensor([coords],
                           requires_grad=False, device=device)
species = consts.species_to_tensor(atoms).to(device).unsqueeze(0)

###############################################################################
# Now let's compute energy and force:
_, energy = model((species, coordinates))
#derivative = torch.autograd.grad(energy.sum(), coordinates)[0]
#force = -derivative

end = timer()

###############################################################################
# And print to see the result:
print(sys.argv[1], -627.509469*(energy.item()), end-start, sep=',' )
