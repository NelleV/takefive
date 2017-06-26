"""Nuclear Landmark model for Plasmodium"""
import os
import sys
import math
import IMP
import IMP.core
import IMP.atom
import IMP.display
import IMP.algebra
import numpy as np
import time
import argparse

state = "RINGS"

parser = argparse.ArgumentParser(description='Run volume exclusion model')
parser.add_argument('seed', metavar='N', type=int,
                    help='seed')
parser.add_argument("--state", default="rings")
args = parser.parse_args()

# Set the seed and the random state
state  = args.state.lower()
random_state = np.random.RandomState(seed=args.seed)

# Nuclear radius - this should change depending on the state of the plasmodium
# Let's work with rings right now
if state == "rings" or state == "r":
    nuclear_rad = 350
elif state in ["schizont", "schizonts", "s"]:
    nuclear_rad = 450
elif state in ["trophozoites", "trophozoite", "t", "troph"]:
    nuclear_rad = 850
else:
    raise ValueError("Unknown state %s" % state)

basedir = os.path.dirname(os.path.realpath(__file__))
outname = os.path.join(basedir, "cluster", str(args.seed),
                       '%s_plasmodium_without_VRSM' % state)
print('output pdb:', outname)
if os.path.exists(outname):
    sys.exit(0)

try:
    os.makedirs(os.path.dirname(outname))
except OSError:
    pass

chr_seq = {"chr1": 640851, "chr2": 947102, "chr3": 1067971, "chr4": 1200490,
           "chr5": 1343557, "chr6": 1418242, "chr7": 1445207, "chr8": 1472805,
           "chr9": 1541735, "chr10": 1687656, "chr11": 2038340,
           "chr12": 2271494,
           "chr13": 2925236, "chr14": 3291936}

# Centromers middle position
chr_cen = {"chr1": 459191, "chr2": 448856, "chr3": 599144, "chr4": 643179,
           "chr5": 456668, "chr6": 479886, "chr7": 810620, "chr8": 300205,
           "chr9": 1243266, "chr10": 936752, "chr11": 833107,
           "chr12": 1283731,
           "chr13": 1169395, "chr14": 1073139}

# Var genes mean position
var_genes = {
    "chr1": [62704, 555068],
    "chr2": [46388, 892057],
    "chr3": [52823, 1020563],
    "chr4": [61257.5, 131090, 159800, 580028,
             961239, 1140018],
    "chr5": [31166, 1332682],
    "chr6": [26775.5,
             732761.5, 1348167.5],
    "chr7": [48681, 561236,  1404322.],
    "chr8": [46573, 448866,  1377043.],
    "chr9": [49477, 1488883.5],
    "chr10": [43735, 1620763.],
    "chr11": [55289, 2021296.],
    "chr12": [36889,
              774639, 1712078.5, 2218861],
    "chr13": [44090, 106156.5,
              2873442],
    "chr14": [17818, 3240144]
    }

chr_pdb = {"chr1": 'c01 A', "chr2": 'c02 B', "chr3": 'c03 C', "chr4": 'c04 D',
           "chr5": 'c05 E', "chr6": 'c06 F', "chr7": 'c07 G', "chr8": 'c08 H',
           "chr9": 'c09 I', "chr10": 'c10 J', "chr11": 'c11 K',
           "chr12": 'c12 L',
           "chr13": 'c13 M', "chr14": 'c14 N'}

chain_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14']


t1 = time.time()
sep = 3200       # 3200 bp separation
chr_bead = {}    # number of beads for each chromosome
nbead = 0
bead_start = {}  # bead label starts of a chr
for i in chr_seq.keys():
    n = chr_seq[i] / sep + 1
    chr_bead[i] = n
    nbead = nbead + n
    bead_start[i] = nbead - n

startChr_pdb = []  # indexing starts from 0
n = 0
startChr_pdb.append(n)
for chrom in chain_list:
    n += chr_bead[chrom]
    startChr_pdb.append(n)


# FUNCTIONS --------------------------------------------------------- start

def bead_id(chr, gpos):
    '''Given chromosome id and genome position, returns bead id'''
    for i in range(chr_bead[chr]):
        if gpos >= i*sep and gpos < (i+1)*sep:
            beadnum = i + bead_start[chr]
            break
    return beadnum


def find_chromosome(bid):
    """ Returns a chromosome id given a bead number"""
    for i in chr_seq.keys():
        if bid < bead_start[i] + chr_bead[i] \
           and bid >= bead_start[i]:
            chrid = i
            break
    return chrid


def find_bead_in_chr(bid):
    """
    Returns a chromosome and bead_order and mid genome position given a
    beadnum
    """
    for i in chr_seq.keys():
        if bid < bead_start[i] + chr_bead[i] \
           and bid >= bead_start[i]:
            order = bid - bead_start[i] + 1  # order starts from 1
            genpos = order * sep - sep / 2
            break
    return i, order, genpos


def pdboutput(name=None):
    X = []
    for i in range(nbead):
        p0 = IMP.core.XYZR(chain.get_particle(i))
        chr_num = filter(lambda k: bead_start[k] <= i, chr_seq.keys())[-1]
        X.append([int(chr_num[3:]), p0.get_x(), p0.get_y(), p0.get_z()])
    return X

def mdstep(model, constraints, t, step):
    sf = IMP.core.RestraintsScoringFunction(constraints)
    o = IMP.atom.MolecularDynamics(model)
    # replace 300 K with 500 K
    #md = IMP.atom.VelocityScalingOptimizerState(xyzr, t, 10)
    md = IMP.atom.VelocityScalingOptimizerState(model, xyzr, t)
    o.set_scoring_function(sf)
    o.add_optimizer_state(md)
    # print 'optimizing with temperature',t,'and',step,'steps'
    s = o.optimize(step)
    o.remove_optimizer_state(md)
    # print 'MD',step,'steps done @',datetime.datetime.now()
    return s


def cgstep(model, constraints, step=1000):
    o = IMP.core.ConjugateGradients(model)
    sf = IMP.core.RestraintsScoringFunction(constraints)
    o.set_scoring_function(sf)
    f = o.optimize(step)
    return f
# FUNCTIONS --------------------------------------------------------- start


# ___________________________ IMP starts _____________________________________
# IMP.set_check_level(IMP.NONE)
IMP.set_log_level(IMP.SILENT)
m = IMP.Model()
r = 15.0
lb = 30.0  # length of bond
kbend = 0.2
contact_dict = {}
xyzr = IMP.core.create_xyzr_particles(m, nbead, r)
chain = IMP.container.ListSingletonContainer(xyzr)

# First beads
corner1 = IMP.algebra.Vector3D(-nuclear_rad, -nuclear_rad, -nuclear_rad)
corner2 = IMP.algebra.Vector3D(nuclear_rad, nuclear_rad, nuclear_rad)
box = IMP.algebra.BoundingBox3D(corner1, corner2)
rdummy = int(random_state.rand() * 10000)
for i in range(rdummy):
    ranvec = IMP.algebra.get_random_vector_in(box)
# ----------------------------------------------------------
# print nbead
for i in range(nbead):
    p0 = chain.get_particle(i)
    IMP.atom.Mass.setup_particle(p0, 1)
    p = IMP.core.XYZR(p0)
    coor = IMP.algebra.get_random_vector_in(box)
    p.set_coordinates(coor)
    p.set_coordinates_are_optimized(True)
    # ch,b,gpos = find_bead_in_chr(i)
    # print ch,b,pdbOrder(ch,gpos)+1,i

# sys.exit()
# ----------------------------------------------------------------------------

# print 'Setting up restraints'
# Create bonds for consecutive beads in a string
bonds = IMP.container.ListSingletonContainer(m)
for id in chr_seq.keys():
    istart = bead_start[id]
    iend = istart + chr_bead[id]
    bp = IMP.atom.Bonded.setup_particle(chain.get_particle(istart))
    for i in range(istart + 1, iend):
        bpr = IMP.atom.Bonded.setup_particle(chain.get_particle(i))
        b = IMP.atom.create_custom_bond(bp, bpr, lb, 2)
        bonds.add(b.get_particle())
        bp = bpr

# Restraint for bonds
bss = IMP.atom.BondSingletonScore(IMP.core.Harmonic(0, 1))
br = IMP.container.SingletonsRestraint(bss, bonds)
#m.add_restraint(br)  # 0

# Set up excluded volume
evr = IMP.core.ExcludedVolumeRestraint(chain)
#m.add_restraint(evr)  # 1


sf = IMP.core.RestraintsScoringFunction([evr, br])
bd = IMP.atom.MolecularDynamics(m)
md = IMP.atom.VelocityScalingOptimizerState(m, xyzr, 10)
o = IMP.core.ConjugateGradients(m)
o.set_scoring_function(sf)

#bd.set_scoring_function(sf)
#bd.add_optimizer_state(md)
#o.set_time_step(10000)
# UPTO here FERHAT ####
# print bead_start

# Set up cap
center = IMP.algebra.Vector3D(0, 0, 0)
ubcell = IMP.core.HarmonicUpperBound(nuclear_rad, 1.0)
sscell = IMP.core.DistanceToSingletonScore(ubcell, center)
rcell = IMP.container.SingletonsRestraint(sscell, chain)
#m.add_restraint(rcell)  # 2

# centromers in radius 300 @-700
centro_rad = 50.0
centro = IMP.algebra.Vector3D(nuclear_rad - centro_rad, 0, 0)
listcentro = IMP.container.ListSingletonContainer(m)
for k in chr_seq.keys():
    j = bead_id(k, chr_cen[k])
    pcen = chain.get_particle(j)
    listcentro.add(pcen)

ubcen = IMP.core.HarmonicUpperBound(centro_rad, 1.0)
sscen = IMP.core.DistanceToSingletonScore(ubcen, centro)
rcentro = IMP.container.SingletonsRestraint(sscen, listcentro)
#m.add_restraint(rcentro)  # 4

constraints = [evr, br, rcell, rcentro]

x_ = pdboutput(outname)
mdstep(m, constraints, 1000000, 500)
mdstep(m, constraints, 500000, 500)
mdstep(m, constraints, 300000, 500)
mdstep(m, constraints, 100000, 500)
mdstep(m, constraints, 5000, 500)
score = cgstep(m, constraints, 1000)
print 'before telo: ', score
                              
#score = cgstep(1000)
#print 'before telo: ', score

# Telomeres near nuclear envelope thickness 50
telo = IMP.container.ListSingletonContainer(m)
# galBead =  bead_id(galpos[0],galpos[1])
# telo.add_particle(chain.get_particle(galBead))
for k in chr_seq.keys():
    j1 = bead_start[k]
    pt = chain.get_particle(j1)
    telo.add(pt)
    j2 = j1 - 1 + chr_bead[k]
    pt = chain.get_particle(j2)
    telo.add(pt)
envelope = nuclear_rad - 50.0
tlb = IMP.core.HarmonicLowerBound(envelope, 1.0)
sst = IMP.core.DistanceToSingletonScore(tlb, center)
rt = IMP.container.SingletonsRestraint(sst, telo)
#m.add_restraint(rt)  # 5

constraints = [evr, br, rcell, rcentro, rt]
mdstep(m, constraints, 500000, 5000)
mdstep(m, constraints, 300000, 5000)
mdstep(m, constraints, 5000, 10000)
score = cgstep(m, constraints, 500)
print 'before angle', score


# outside centro sphere
# lbcen = IMP.coreHarmonicLowerBound(centro_rad,1.0)
# sstc = IMP.core.DistanceToSingletonScore(lbcen,centro)
# rtc = IMP.container.SingletonsRestraint(sstc,telo)
# m.add_restraint(rtc) #6

# Set up Nucleolis #
# -------------------

print 'High temp MD in nuc ...'


# Angle Restraint
angle = math.pi
angle_set = []
noangle = [i for i in bead_start.values()]  # do not apply angle restraints
ars = []
for i in range(nbead - 1):
    ieval = i + 1
    if ieval in noangle:
        continue
    elif i in noangle:
        continue
    else:
        d1 = chain.get_particle(i - 1)
        d2 = chain.get_particle(i)
        d3 = chain.get_particle(i + 1)
        pot = IMP.core.Harmonic(angle, kbend)
        ar = IMP.core.AngleRestraint(m, pot, d1, d2, d3)
        ars.append(ar)
        angle_set.append(ar)


constraints = [evr, br, rcell, rcentro, rt] + ars

mdstep(m, constraints, 500000, 5000)
mdstep(m, constraints, 300000, 5000)
mdstep(m, constraints, 5000, 10000)
score = cgstep(m, constraints, 500)



X = pdboutput(outname)

# -------------------------

# mdstep(1000,1000)
# score=cgstep(1000)
# print '\nFinal score with remove angle restraint is: ',score

# name='final_without_angle'
# output(chain,nbead,name)
t2 = time.time()

print 'time spend is ', t2 - t1, ' s'
