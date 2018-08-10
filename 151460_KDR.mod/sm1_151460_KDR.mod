NEURON
{
  SUFFIX KDR 
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

  an = 0.09719524079113401     (/mV) 
  bn = -1.6609239195555587     (1) 
  vhn = -34.982592350795635     (mV) 
  An = 10.304358320846472     (/ms) 
  b1n = -0.05310473475810375     (/mV) 
  c1n = 0.00042036484799905374     (/mV2) 
  d1n = -1.3566286761561015e-06     (/mV3) 
  b2n = -0.06571737982661532     (/mV) 
  c2n = -0.0009693772810752884     (/mV2) 
  d2n = -7.523742941388525e-06     (/mV3) 
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