NEURON
{
  SUFFIX kpkj2 
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

  an = 0.04901943209193984     (/mV) 
  bn = -1.7156767419491925     (1) 
  vhn = -51.96775337204331     (mV) 
  An = 10.411259295168522     (/ms) 
  b1n = -0.09973789514940559     (/mV) 
  c1n = 0.0011005751325551727     (/mV2) 
  d1n = -4.290558721139837e-06     (/mV3) 
  b2n = -0.08685124622401144     (/mV) 
  c2n = -0.002202153601083006     (/mV2) 
  d2n = -2.3422290346668373e-05     (/mV3) 
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