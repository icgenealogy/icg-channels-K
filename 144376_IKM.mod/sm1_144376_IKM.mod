NEURON
{
  SUFFIX im 
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

  am = 0.099967538074007     (/mV) 
  bm = -3.4988379641369067     (1) 
  vhm = -46.967104130922955     (mV) 
  Am = 506.49974645245345     (/ms) 
  b1m = -0.04993641333916172     (/mV) 
  c1m = -1.0467574899838353e-06     (/mV2) 
  d1m = 6.8079887854237676e-09     (/mV3) 
  b2m = -0.050075860598188214     (/mV) 
  c2m = -1.6545016820869189e-06     (/mV2) 
  d2m = -1.3173370950688902e-08     (/mV3) 
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