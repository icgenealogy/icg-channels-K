NEURON
{
  SUFFIX Kirm 
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

  am = -0.09999989455104703     (/mV) 
  bm = 9.999993136233048     (1) 
  vhm = -320.9799344550441     (mV) 
  Am = 9.246595584919156     (/ms) 
  b1m = 4.260296867619112     (/mV) 
  c1m = -0.021316494229465827     (/mV2) 
  d1m = 2.656334251403521e-05     (/mV3) 
  b2m = 0.05501095435178393     (/mV) 
  c2m = -0.00014659947754466676     (/mV2) 
  d2m = 1.1906553775227271e-07     (/mV3) 
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