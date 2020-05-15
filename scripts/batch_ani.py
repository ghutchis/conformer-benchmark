# -*- coding: utf-8 -*-
"""
Computing energies using TorchANI
"""

from __future__ import print_function

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

import torch
import torchani

from timeit import default_timer as timer

import sys
import os

# Manually specify the device we want TorchANI to run:
# benchmarks were run on a server without CUDA so CPU is the real default
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#######
# We're going to read in an XYZ file specified on the command-line
coord_sets = []
if len(sys.argv) > 2:
    for name in sys.argv[2:]:
        with open(name) as xyz_file:
            atoms = []
            coords = []
            num_atoms = int(xyz_file.readline().rstrip())
            xyz_file.readline() # skip the smiles or title line
            for i in range(num_atoms):
                element, x, y, z = xyz_file.readline().split()
                atoms.append(element)
                coords.append([float(x),float(y),float(z)])
            coord_sets.append(coords)

###############################################################################
# Let's now load the built-in ANI-1 models. The builtin ANI-1ccx contains 8
# models trained with diffrent initialization. Predicting the energy and force
# using the average of the 8 models outperform using a single model, so it is
# always recommended to use an ensemble, unless the speed of computation is an
# issue in your application.

# (not sure if these need the level shifters...)
ani1x = torchani.models.ANI1x()
ani1ccx = torchani.models.ANI1ccx()

# path to ANI-2x subdirectories
path = '/Users/ghutchis/Devel/torchani/torchani/resources/ani-2x'

# this whole pile is to load ANI-2x (since it's not public)
#const_file = path + '/rHCNOSFCl-5.1R_16-3.5A_a8-4.params'  # noqa: E501
#consts = torchani.neurochem.Constants(const_file)
#aev_computer = torchani.AEVComputer(**consts)
#sae_file = path + '/sae_linfit.dat'  # noqa: E501
#energy_shifter = torchani.neurochem.load_sae(sae_file)
#model_prefix = path + '/train'  # noqa: E501
#ensemble = torchani.neurochem.load_model_ensemble(consts.species, model_prefix, 8)  # noqa: E501

#ani2x = torch.nn.Sequential(aev_computer, ensemble, energy_shifter)
# done (that's ANI-2x)

model = ani1ccx

###############################################################################
# Now let's define the coordinate and species. If you just want to compute the
# energy and force for a single structure like in this example, you need to
# make the coordinate tensor has shape ``(1, Na, 3)`` and species has shape
# ``(1, Na)``, where ``Na`` is the number of atoms in the molecule, the
# preceding ``1`` in the shape is here to support batch processing like in
# training. If you have ``N`` different structures to compute, then make it
# ``N``.
species = model.species_to_tensor(atoms).to(device).unsqueeze(0)

# Now let's compute energy and force:
for coords in coord_sets:
    start = timer()
    coordinates = torch.tensor([coords], requires_grad=False, device=device)
    _, energy = model((species, coordinates))
    end = timer()

    # Calculating an atomization energy from ANI was essentially impossible
    # We're converting to kcal/mol and multiplying by -1.0 to get the right sign
    print(sys.argv[1], -627.509469*(energy.item()), end-start, sep=',' )
