NEURON
{
  SUFFIX KAfm 
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

  ah = -0.13157856874027052     (/mV) 
  bh = 9.263133713735186     (1) 
  vhh = -68.34435981863342     (mV) 
  Ah = 12.660000007925616     (/ms) 
  b1h = -3.4961349863121866e-09     (/mV) 
  c1h = -2.0648845158487e-09     (/mV2) 
  d1h = -3.269319039542681e-11     (/mV3) 
  b2h = -3.5203637451291676e-09     (/mV) 
  c2h = -2.062626150278558e-09     (/mV2) 
  d2h = -3.271447379471349e-11     (/mV3) 

  am = 0.1333331054859963     (/mV) 
  bm = -4.413319532547033     (1) 
  vhm = -61.286462867376045     (mV) 
  Am = 0.5200000003979275     (/ms) 
  b1m = -2.7276648320281297e-09     (/mV) 
  c1m = 5.546166329054749e-10     (/mV2) 
  d1m = 4.5278765934308743e-11     (/mV3) 
  b2m = -2.731955978592326e-09     (/mV) 
  c2m = 5.566312776564243e-10     (/mV2) 
  d2m = 4.5257903285997825e-11     (/mV3) 
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
  g = gbar*h*m
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