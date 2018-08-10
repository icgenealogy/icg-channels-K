NEURON
{
  SUFFIX kdrcurrent 
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

  an = 0.10999480678833093     (/mV) 
  bn = 1.4300592057176864     (1) 
  vhn = 20.607103908380267     (mV) 
  An = 38.07299499795778     (/ms) 
  b1n = -0.10327700905145538     (/mV) 
  c1n = 0.00046254016332927686     (/mV2) 
  d1n = 8.044666244423435e-06     (/mV3) 
  b2n = -0.005659597528581894     (/mV) 
  c2n = 0.00037782775606807486     (/mV2) 
  d2n = 1.8465417479665448e-06     (/mV3) 
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
  g = gbar
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