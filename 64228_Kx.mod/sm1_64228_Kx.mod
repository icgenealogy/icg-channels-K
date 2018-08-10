NEURON
{
  SUFFIX Kx 
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

  anKx = 0.17559567135846144     (/mV) 
  bnKx = -8.763973707470823     (1) 
  vhnKx = -49.90677146903192     (mV) 
  AnKx = 1515.1543330391435     (/ms) 
  b1nKx = 0.0877698626882023     (/mV) 
  c1nKx = 1.5404028536385472e-06     (/mV2) 
  d1nKx = 1.4442821913797432e-08     (/mV3) 
  b2nKx = 0.08766521725462861     (/mV) 
  c2nKx = 1.6892249048056357e-06     (/mV2) 
  d2nKx = -1.9286398812434923e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nKxInf 
  nKxTau 
}

STATE
{
  nKx
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*nKx
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  nKx' = (nKxInf - nKx) / nKxTau 
}

INITIAL
{
  rates(v)
  nKx = nKxInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nKxInf = 1/(1 + exp(-anKx*v + bnKx)) 
    nKxTau = AnKx / ( exp(-(b1nKx*(v-vhnKx) + c1nKx*(v-vhnKx)^2 + d1nKx*(v-vhnKx)^3)) + exp((b2nKx*(v-vhnKx) + c2nKx*(v-vhnKx)^2 + d2nKx*(v-vhnKx)^3)) ) 


  UNITSON
}