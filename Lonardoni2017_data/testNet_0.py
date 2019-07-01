__author__ = 'dlonardoni'

import matplotlib as m
m.use('Qt5Agg')

import pylab as plt
import cellCultureNet as Net

plt.ion()
Net.StartSmallNetTest()
plt.show()