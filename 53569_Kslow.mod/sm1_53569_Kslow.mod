NEURON
{
  SUFFIX Ks 
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

  an = 0.09068699167586286     (/mV) 
  bn = -2.611502116354381     (1) 
  vhn = -51.33914970431323     (mV) 
  An = 28.694354212103793     (/ms) 
  b1n = -0.04160251969843303     (/mV) 
  c1n = 0.0002662222923843507     (/mV2) 
  d1n = -7.2001513918378e-07     (/mV3) 
  b2n = -0.06958104337497674     (/mV) 
  c2n = -0.0011891540763519215     (/mV2) 
  d2n = -1.1130659897753356e-05     (/mV3) 
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