NEURON
{
  SUFFIX km 
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

  an = 0.11105427483423269     (/mV) 
  bn = -3.3315909210568617     (1) 
  vhn = -34.86901953639041     (mV) 
  An = 33.82909619496914     (/ms) 
  b1n = -0.04508861048879203     (/mV) 
  c1n = 0.0003042542228761562     (/mV2) 
  d1n = -8.619027078581126e-07     (/mV3) 
  b2n = -0.06381649957169595     (/mV) 
  c2n = -0.0007867578962415032     (/mV2) 
  d2n = -4.4063306653217804e-06     (/mV3) 
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