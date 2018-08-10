NEURON
{
  SUFFIX ka 
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

  aa = 0.04999993416493523     (/mV) 
  ba = -1.4999953428328596     (1) 
  vha = -229.3119005287641     (mV) 
  Aa = 0.509995008551842     (/ms) 
  b1a = 1.7939121066805808e-07     (/mV) 
  c1a = -9.188827806480189e-10     (/mV2) 
  d1a = 1.424377435865485e-12     (/mV3) 
  b2a = -25.86100719321309     (/mV) 
  c2a = 0.3192706783443542     (/mV2) 
  d2a = -0.0009860870400413906     (/mV3) 

  ab = -0.16666458508998017     (/mV) 
  bb = 13.33318814134198     (1) 
  vhb = -249.27712281119668     (mV) 
  Ab = 33.8962993114201     (/ms) 
  b1b = -0.005211796328878517     (/mV) 
  c1b = 1.2123195898987042e-05     (/mV2) 
  d1b = -9.910099960721682e-09     (/mV3) 
  b2b = -0.00023248592574012352     (/mV) 
  c2b = -3.2595026585493385e-05     (/mV2) 
  d2b = 4.388795819388664e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  aInf 
  aTau 
  bInf 
  bTau 
}

STATE
{
  a
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*a*a*a*b
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  a' = (aInf - a) / aTau 
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  a = aInf 
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    aInf = 1/(1 + exp(-aa*v + ba)) 
    aTau = Aa / ( exp(-(b1a*(v-vha) + c1a*(v-vha)^2 + d1a*(v-vha)^3)) + exp((b2a*(v-vha) + c2a*(v-vha)^2 + d2a*(v-vha)^3)) ) 

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}