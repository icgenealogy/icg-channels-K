NEURON
{
  SUFFIX emdk 
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

  an = 0.09090601241811166     (/mV) 
  bn = 1.2727369165751     (1) 
  vhn = 26.6043038134597     (mV) 
  An = 47.83142504799483     (/ms) 
  b1n = 0.03111539688004845     (/mV) 
  c1n = -1.3119269719336522e-06     (/mV2) 
  d1n = -6.115672697319347e-09     (/mV3) 
  b2n = 0.0335213689046012     (/mV) 
  c2n = -2.7671765858667993e-06     (/mV2) 
  d2n = 1.135073932181583e-08     (/mV3) 
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