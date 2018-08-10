NEURON
{
  SUFFIX Kdr 
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

  an = 0.046834510013836504     (/mV) 
  bn = -0.9012470838569978     (1) 
  vhn = 6.019553175384406     (mV) 
  An = 3.93435818765325     (/ms) 
  b1n = -0.020397117079041384     (/mV) 
  c1n = -2.9349478522508666e-05     (/mV2) 
  d1n = 1.1185079037284351e-06     (/mV3) 
  b2n = -0.003990640619077527     (/mV) 
  c2n = 9.703459473185021e-05     (/mV2) 
  d2n = 6.588074716038998e-07     (/mV3) 
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