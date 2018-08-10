NEURON
{
  SUFFIX iA 
  USEION k READ ek WRITE ik 
  RANGE gbar, g, ik
  GLOBAL ek
}

UNITS
{
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)
}

PARAMETER
{
  gbar = 1 (S/cm2)

  an = 0.03448271053806346     (/mV) 
  bn = -1.6206865226498475     (1) 
  vhn = 17.651949278429466     (mV) 
  An = 0.2376667052071567     (/ms) 
  b1n = -0.010764444557234352     (/mV) 
  c1n = -3.8900862643213904e-05     (/mV2) 
  d1n = 3.2189344277843024e-07     (/mV3) 
  b2n = 0.010760583957543621     (/mV) 
  c2n = -0.0003683472982187285     (/mV2) 
  d2n = -1.2969988160184496e-06     (/mV3) 

  al = -0.0999990625648875     (/mV) 
  bl = 6.599953501326262     (1) 
  vhl = -70.74568503186194     (mV) 
  Al = 30.176300553227218     (/ms) 
  b1l = -0.08866574944720741     (/mV) 
  c1l = 0.0008811314729856003     (/mV2) 
  d1l = -2.6974036529181667e-06     (/mV3) 
  b2l = -0.07750989744815891     (/mV) 
  c2l = -0.002025252929510281     (/mV2) 
  d2l = -4.867879399577431e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
  lInf 
  lTau 
}

STATE
{
  n
  l
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*l
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
  l' = (lInf - l) / lTau 
}

INITIAL
{
  rates(v)
  n = nInf 
  l = lInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 


  UNITSON
}