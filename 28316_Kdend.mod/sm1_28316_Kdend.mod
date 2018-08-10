NEURON
{
  SUFFIX Kdend 
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

  an = 0.05184634221022714     (/mV) 
  bn = -0.6935143274206388     (1) 
  vhn = -11.258152422890369     (mV) 
  An = 1.5179380812691226     (/ms) 
  b1n = -0.01908919689392665     (/mV) 
  c1n = -0.0001059442711218748     (/mV2) 
  d1n = 1.1977040623662027e-06     (/mV3) 
  b2n = -0.009443512275483713     (/mV) 
  c2n = 6.684988615339075e-05     (/mV2) 
  d2n = 9.61722318440557e-07     (/mV3) 
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
  g = gbar*n*n*n*n
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