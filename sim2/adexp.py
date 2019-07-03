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
def createcell(section, **kwargs):
    cell = h.AdExp(0.5, sec=section) # Create a new Izhikevich neuron at location 0 (doesn't matter where) in h.Section() "section"
    cell.rpeso = 30
    cell.mNMDA = 0*1
    cell.v0_block = -55
    cell.k_block = 5
    cell.mAMPA = 1
    cell.iEXT = 50
    cell.tau_w = 280
    cell.G_l = 10
    cell.a = 2
    cell.C = 281
    cell.E_l = -70.6
    cell.V_thre = -50.4
    cell.Delta_T = 2
    cell.V_reset = -70.6
    cell.b = 40
    return cell
