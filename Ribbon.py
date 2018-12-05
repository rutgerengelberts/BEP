import numpy as np
import kwant
from matplotlib import pyplot
from math import atan2, pi, sqrt

pot = 0.1;
W = 5;
L = 5.5;

graphene = kwant.lattice.honeycomb()

def rect(pos):
    x, y = pos
    return -W < x < W and -L < y < L

def potential(site):
    (x, y) = site.pos
    d = y * cos_30 + x * sin_30
    return pot * tanh(d / w)

syst = kwant.Builder()
syst[graphene.shape(rect, (0, 0))] = potential
syst[graphene.neighbors()] = -1

def lead(pos):
    x, y = pos
    return -L < y < L

sym0 = kwant.TranslationalSymmetry(graphene.vec((-1, 0)))
leadleft = kwant.Builder(sym0)
leadleft[graphene.shape(lead, (0, 0))] = -pot
leadleft[graphene.neighbors()] = -1

sym1 = kwant.TranslationalSymmetry(graphene.vec((1, 0)))
leadright = kwant.Builder(sym1)
leadright[graphene.shape(lead, (0, 0))] = -pot
leadright[graphene.neighbors()] = -1

syst.attach_lead(leadleft)
syst.attach_lead(leadright)

kwant.plot(syst)


syst = syst.finalized()
momenta = [-2*pi + 0.02 * 2*pi * i for i in range(101)]
bands = kwant.physics.Bands(syst.leads[0], momenta)
energies = [bands(k) for k in momenta]
pyplot.figure()
pyplot.plot(momenta, energies)
pyplot.xlabel("Momentum [(lattice constant)^-1]")
pyplot.ylabel("Energy [t]")
pyplot.show()
















