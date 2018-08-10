NEURON
{
  SUFFIX SKv3_1 
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

  am = 0.10309259748818682     (/mV) 
  bm = 1.9278445847207941     (1) 
  vhm = -39.412963127089014     (mV) 
  Am = 4.332983764497996     (/ms) 
  b1m = -0.0012148404035759715     (/mV) 
  c1m = 7.474636284453255e-06     (/mV2) 
  d1m = -1.7285541363474698e-08     (/mV3) 
  b2m = -0.022015098296621335     (/mV) 
  c2m = -1.4693177681754831e-08     (/mV2) 
  d2m = 4.541537847942305e-09     (/mV3) 
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