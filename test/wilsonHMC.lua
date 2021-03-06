--doTrace = true
nf = 2
nx = 4
nt = 8
seed = 3190
beta = 5.7
gact = { type="symanzik_1loop", u0=0.8695 }
gintalg = { type="leapfrog" }
fintalg = { type="leapfrog" }
--hmcmasses = { 0, 0.1 }
--hmcmasses = { 0.38596491228070175438 }
hmcmasses = { 0 }
--rho = 0.2
tau = 0.1
ngsteps = 20
nfsteps = { 10 }
faresid = 1e-8
mdresid = 1e-8
warmup = 1
ntraj = 1
nlats = 2

require 'wilson'

TESTPAT("^plaq")
TESTPAT("^S")
TESTPAT("deltaS")
TESTPAT("^1")
TESTPAT("MEAS")
