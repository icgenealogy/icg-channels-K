NEURON
{
  SUFFIX Ikt1m4h 
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

  ah = -0.2717410750893756     (/mV) 
  bh = 10.752827853352109     (1) 
  vhh = -38.72825618461124     (mV) 
  Ah = 87.38157195991919     (/ms) 
  b1h = -0.12637081915811982     (/mV) 
  c1h = 0.0016274177874946677     (/mV2) 
  d1h = -6.376217860237394e-06     (/mV3) 
  b2h = -0.3682732847218359     (/mV) 
  c2h = -0.01184137228292776     (/mV2) 
  d2h = -0.0001103491444525492     (/mV3) 

  am = 0.05387923309410741     (/mV) 
  bm = -2.106675664575642     (1) 
  vhm = -32.62095594346716     (mV) 
  Am = 1.391290939316626     (/ms) 
  b1m = -0.0707259784105378     (/mV) 
  c1m = 0.0008688539149073047     (/mV2) 
  d1m = -3.3275172010706052e-06     (/mV3) 
  b2m = -0.08972721702994081     (/mV) 
  c2m = -0.0017357996486448523     (/mV2) 
  d2m = -1.1822107750388201e-05     (/mV3) 
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
  g = gbar*h*m*m*m*m
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