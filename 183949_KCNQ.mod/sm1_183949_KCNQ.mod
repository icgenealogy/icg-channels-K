NEURON
{
  SUFFIX KCNQ 
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

  am = 0.0512831095408758     (/mV) 
  bm = -3.1283545046788563     (1) 
  vhm = -44.6524701060197     (mV) 
  Am = 86.65832004683554     (/ms) 
  b1m = -0.04975584069665789     (/mV) 
  c1m = 0.00032798313201176687     (/mV2) 
  d1m = -6.328299708037163e-07     (/mV3) 
  b2m = -0.010334154252062284     (/mV) 
  c2m = 0.0001105917066282893     (/mV2) 
  d2m = 1.2084245934944867e-07     (/mV3) 
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