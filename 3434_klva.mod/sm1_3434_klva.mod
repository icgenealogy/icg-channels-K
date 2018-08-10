NEURON
{
  SUFFIX klva 
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

  an = 0.11729684631734283     (/mV) 
  bn = -7.200288361493359     (1) 
  vhn = -57.641338164337846     (mV) 
  An = 2.0771873201251427     (/ms) 
  b1n = -0.05367429404911354     (/mV) 
  c1n = 0.00011189966018053846     (/mV2) 
  d1n = -2.505805372124558e-07     (/mV3) 
  b2n = -0.05441987250808852     (/mV) 
  c2n = 0.0004135438300040574     (/mV2) 
  d2n = 3.7462344799163923e-06     (/mV3) 
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