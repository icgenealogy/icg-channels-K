NEURON
{
  SUFFIX kdrJonas 
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

  am = 0.056254590949879404     (/mV) 
  bm = -2.8702504185428057     (1) 
  vhm = -73.86025321886802     (mV) 
  Am = 11.566691077442446     (/ms) 
  b1m = -0.03731813872782305     (/mV) 
  c1m = 0.00015560377125469298     (/mV2) 
  d1m = -9.5457316832068e-08     (/mV3) 
  b2m = -0.0317023720116892     (/mV) 
  c2m = -0.0001779092104991066     (/mV2) 
  d2m = 2.1628213206443936e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}