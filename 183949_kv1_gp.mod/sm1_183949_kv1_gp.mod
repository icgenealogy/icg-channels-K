NEURON
{
  SUFFIX kv1_gp 
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

  ah = -0.04604480642918483     (/mV) 
  bh = 1.5298259728787709     (1) 
  vhh = -75.80410052463957     (mV) 
  Ah = 395.901476519724     (/ms) 
  b1h = 0.0020293775932114086     (/mV) 
  c1h = 0.00014136528405862218     (/mV2) 
  d1h = -8.453696136987559e-07     (/mV3) 
  b2h = -0.0020289034484901633     (/mV) 
  c2h = 0.00013035437060827308     (/mV2) 
  d2h = -5.759928428337937e-07     (/mV3) 

  am = 0.06248412333735882     (/mV) 
  bm = -1.6865547423769123     (1) 
  vhm = -17.64020878102342     (mV) 
  Am = 26.574492269808733     (/ms) 
  b1m = -0.07351041325088768     (/mV) 
  c1m = -0.0003128693544295813     (/mV2) 
  d1m = 8.301861312867361e-06     (/mV3) 
  b2m = 0.0002906100463122681     (/mV) 
  c2m = 0.0005024620728429028     (/mV2) 
  d2m = -1.8378046511891378e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}