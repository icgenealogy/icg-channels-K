NEURON
{
  SUFFIX Khcvode 
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

  am = -0.14284760309695477     (/mV) 
  bm = 11.142221759670113     (1) 
  vhm = -32.65823390129214     (mV) 
  Am = 92.26000001099862     (/ms) 
  b1m = 3.0942528479985416e-09     (/mV) 
  c1m = 1.004768719706702e-09     (/mV2) 
  d1m = 6.15433784908847e-11     (/mV3) 
  b2m = 3.1441440896244175e-09     (/mV) 
  c2m = 1.0050117039434909e-09     (/mV2) 
  d2m = 6.152920743368169e-11     (/mV3) 
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