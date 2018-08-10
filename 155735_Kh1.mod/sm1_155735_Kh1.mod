NEURON
{
  SUFFIX Kh1 
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

  am = -0.10587071070564114     (/mV) 
  bm = 9.200787492183148     (1) 
  vhm = -75.38429253312358     (mV) 
  Am = 15.212456899862357     (/ms) 
  b1m = 0.0013792606482909392     (/mV) 
  c1m = 3.4031496826073642e-06     (/mV2) 
  d1m = 1.1249486955897718e-08     (/mV3) 
  b2m = 0.0013177503342883003     (/mV) 
  c2m = 2.71876047484275e-06     (/mV2) 
  d2m = -4.976312026628765e-09     (/mV3) 
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