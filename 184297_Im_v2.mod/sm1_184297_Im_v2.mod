NEURON
{
  SUFFIX Im_v2 
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

  am = 0.2297239551211746     (/mV) 
  bm = -11.026840840025846     (1) 
  vhm = -34.780239931969675     (mV) 
  Am = 66.29181247338572     (/ms) 
  b1m = -0.11726393940245454     (/mV) 
  c1m = 0.0020902373428630322     (/mV2) 
  d1m = -9.623683851177296e-06     (/mV3) 
  b2m = -0.013530855346992414     (/mV) 
  c2m = 0.0007335854375267874     (/mV2) 
  d2m = -4.108976421864169e-06     (/mV3) 
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