NEURON
{
  SUFFIX Kdrax 
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

  an = 0.09171100305975641     (/mV) 
  bn = -3.4026118656548396     (1) 
  vhn = -61.26047992470738     (mV) 
  An = 3.4636926670776274     (/ms) 
  b1n = -0.037582271954105455     (/mV) 
  c1n = 0.0001773240760039677     (/mV2) 
  d1n = -2.0311602496111009e-07     (/mV3) 
  b2n = -0.07179317472025323     (/mV) 
  c2n = -0.0013963595346566252     (/mV2) 
  d2n = -1.551991826599395e-05     (/mV3) 
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