NEURON
{
  SUFFIX kd 
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

  an = -0.11230885024266093     (/mV) 
  bn = 3.7027776952192357     (1) 
  vhn = -36.35389243863696     (mV) 
  An = 91.75088685120075     (/ms) 
  b1n = 0.0884135898835971     (/mV) 
  c1n = 3.106989584357738e-05     (/mV2) 
  d1n = -5.0398017311917564e-06     (/mV3) 
  b2n = 0.022543801694376357     (/mV) 
  c2n = 0.00019535660864015134     (/mV2) 
  d2n = -1.0645710630651833e-06     (/mV3) 
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