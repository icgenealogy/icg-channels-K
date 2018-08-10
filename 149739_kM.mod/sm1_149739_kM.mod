NEURON
{
  SUFFIX kM 
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

  ax = 0.20005452428432974     (/mV) 
  bx = -7.002266719569006     (1) 
  vhx = -50.93088049541514     (mV) 
  Ax = 451.06925209122966     (/ms) 
  b1x = -0.024970721441049393     (/mV) 
  c1x = -3.900743467359412e-07     (/mV2) 
  d1x = 1.7346010928187442e-09     (/mV3) 
  b2x = -0.04999684938008936     (/mV) 
  c2x = 6.269880144404626e-07     (/mV2) 
  d2x = 1.0814597059808854e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  xInf 
  xTau 
}

STATE
{
  x
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*x
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  x' = (xInf - x) / xTau 
}

INITIAL
{
  rates(v)
  x = xInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    xInf = 1/(1 + exp(-ax*v + bx)) 
    xTau = Ax / ( exp(-(b1x*(v-vhx) + c1x*(v-vhx)^2 + d1x*(v-vhx)^3)) + exp((b2x*(v-vhx) + c2x*(v-vhx)^2 + d2x*(v-vhx)^3)) ) 


  UNITSON
}