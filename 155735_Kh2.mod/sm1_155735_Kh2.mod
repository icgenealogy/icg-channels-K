NEURON
{
  SUFFIX Kh2 
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

  am = -0.060334158474779     (/mV) 
  bm = 7.345774765565939     (1) 
  vhm = -28.156811267270122     (mV) 
  Am = 73.62043549990936     (/ms) 
  b1m = -0.01077212720651434     (/mV) 
  c1m = 5.881948580456654e-05     (/mV2) 
  d1m = -1.1307449635965494e-07     (/mV3) 
  b2m = -0.010768683005389308     (/mV) 
  c2m = -5.719423981852185e-05     (/mV2) 
  d2m = -9.936124602306829e-08     (/mV3) 
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