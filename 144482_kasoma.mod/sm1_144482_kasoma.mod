NEURON
{
  SUFFIX kasoma 
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

  ah = -0.16666523753663381     (/mV) 
  bh = 12.99990047523851     (1) 
  vhh = -230.1983040519814     (mV) 
  Ah = 21.963314619987628     (/ms) 
  b1h = 0.00026247859202995074     (/mV) 
  c1h = 4.105207108203303e-05     (/mV2) 
  d1h = -3.4223780318475744e-08     (/mV3) 
  b2h = 0.006440899911561653     (/mV) 
  c2h = -1.7216923422783574e-05     (/mV2) 
  d2h = 1.5884578726662982e-08     (/mV3) 

  am = 0.1176464778188413     (/mV) 
  bm = -7.058778104944371     (1) 
  vhm = -60.67643255034253     (mV) 
  Am = 2.3549870641179678     (/ms) 
  b1m = -0.059671819446983625     (/mV) 
  c1m = 0.0004732710526316983     (/mV2) 
  d1m = -1.2472265235636963e-06     (/mV3) 
  b2m = -0.07195941647852114     (/mV) 
  c2m = -0.0003250803030413409     (/mV2) 
  d2m = 3.876960289219264e-06     (/mV3) 
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