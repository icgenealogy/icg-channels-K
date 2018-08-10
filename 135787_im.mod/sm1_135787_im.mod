NEURON
{
  SUFFIX iM 
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

  am = 0.10000204005125472     (/mV) 
  bm = -3.500121053700705     (1) 
  vhm = -35.06596401983419     (mV) 
  Am = 68.7764145591819     (/ms) 
  b1m = 0.05006589464017311     (/mV) 
  c1m = -1.428002825974796e-06     (/mV2) 
  d1m = -3.353439715621131e-08     (/mV3) 
  b2m = 0.04974172483763199     (/mV) 
  c2m = 6.915405667389778e-06     (/mV2) 
  d2m = -6.246341887922063e-08     (/mV3) 
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
  g = gbar*m
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