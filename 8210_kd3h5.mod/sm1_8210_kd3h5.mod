NEURON
{
  SUFFIX kd3 
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

  an = 0.11104961251521618     (/mV) 
  bn = -0.08006039935375078     (1) 
  vhn = 3.2888825482716215     (mV) 
  An = 1.9712038879616867     (/ms) 
  b1n = -0.0823736326233949     (/mV) 
  c1n = 0.0007518950544753914     (/mV2) 
  d1n = -2.6767768354631e-06     (/mV3) 
  b2n = -0.025659434122993135     (/mV) 
  c2n = -0.0001441635916213034     (/mV2) 
  d2n = -3.9261676527100766e-07     (/mV3) 
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