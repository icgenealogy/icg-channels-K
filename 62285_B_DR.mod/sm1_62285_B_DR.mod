NEURON
{
  SUFFIX B_DR 
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

  an = 0.08256607284337843     (/mV) 
  bn = -2.6286979930804955     (1) 
  vhn = -40.03662303128411     (mV) 
  An = 0.5534073920462419     (/ms) 
  b1n = 0.049430933369866197     (/mV) 
  c1n = 0.0003830543887837087     (/mV2) 
  d1n = 2.3093862256112937e-06     (/mV3) 
  b2n = 0.046981834690423004     (/mV) 
  c2n = -0.0003509463729685873     (/mV2) 
  d2n = 1.0926447480764486e-06     (/mV3) 
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