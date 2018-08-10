NEURON
{
  SUFFIX K_T 
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

  ah = -0.09999906256488898     (/mV) 
  bh = 6.599953501326365     (1) 
  vhh = -70.74568503186194     (mV) 
  Ah = 30.176300553227218     (/ms) 
  b1h = -0.08866574944720741     (/mV) 
  c1h = 0.0008811314729856003     (/mV2) 
  d1h = -2.6974036529181667e-06     (/mV3) 
  b2h = -0.07750989744815891     (/mV) 
  c2h = -0.002025252929510281     (/mV2) 
  d2h = -4.867879399577431e-05     (/mV3) 

  am = 0.034482710538063926     (/mV) 
  bm = -1.6206865226498677     (1) 
  vhm = 17.651949278429466     (mV) 
  Am = 0.2376667052071567     (/ms) 
  b1m = -0.010764444557234352     (/mV) 
  c1m = -3.8900862643213904e-05     (/mV2) 
  d1m = 3.2189344277843024e-07     (/mV3) 
  b2m = 0.010760583957543621     (/mV) 
  c2m = -0.0003683472982187285     (/mV2) 
  d2m = -1.2969988160184496e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m*m*m*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}