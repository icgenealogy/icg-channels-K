NEURON
{
  SUFFIX Potassium 
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

  an = 0.0884922278926221     (/mV) 
  bn = -1.5308412370266447     (1) 
  vhn = -56.62693455784713     (mV) 
  An = 6.0005078376623695     (/ms) 
  b1n = -0.00018466029524042842     (/mV) 
  c1n = 2.2647950913083957e-05     (/mV2) 
  d1n = -1.3040239857875328e-07     (/mV3) 
  b2n = -9.51252861532245e-05     (/mV) 
  c2n = 1.973862465841438e-05     (/mV2) 
  d2n = -1.1548906537838679e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
}

STATE
{
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}