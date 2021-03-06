import pstats, cProfile

import pyximport
pyximport.install()

import MutualInductance

cProfile.runctx("MutualInductance.mutualInductance()", globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()
