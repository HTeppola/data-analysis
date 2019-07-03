"""
ADEXP

Python wrapper for an AdExp neuron.

Usage example:
    from neuron import h
    from adexp import createcell
    dummy = h.Section()
    cell = createcell(dummy)

Version: 2019jul03 by cliffk
"""

from neuron import h # Open NEURON

## Create basic AdExp neuron with default parameters
def createcell(section, cellid, **kwargs):
    cell = h.AdExp(0.5, sec=section) # Create a new Izhikevich neuron at location 0 (doesn't matter where) in h.Section() "section"
#    cell.rpeso = rndParam['peso']
#    cell.mNMDA = modPar['mNMDA']
#    cell.v0_block = modPar['vblock']
#    cell.k_block = modPar['kb']
#    cell.mAMPA = modPar['mAMPA']
#    cell.iEXT = modPar['iEXT']
#    cell.tau_w = modPar['tau_w'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
#    cell.G_l = modPar['G_l'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
#    cell.a = modPar['a'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
#    cell.C = 281
#    cell.E_l = modPar['E_l'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
#    cell.V_thre = -50.4
#    cell.Delta_T = 2
#    cell.V_reset = -70.6
#    cell.b = modPar['b'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
#    cell.cellid = cellid # Cell ID for keeping track which cell this is
    return cell

#modParexc = {'tau_w': 280, 'G_l': 10, 'a': 2, 'b': 40, 'C': , 'E_l': -70.6, 'Delta_T': 2,
#             'V_reset': , 'type': 'undef', 'iEXT': 50, 'amp': 0}