NEURON
{
  SUFFIX Kdrsd 
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

  an = 0.09690275663581276     (/mV) 
  bn = -2.0769606292157565     (1) 
  vhn = -40.4870522416352     (mV) 
  An = 8.143361808276003     (/ms) 
  b1n = 0.06705755841381328     (/mV) 
  c1n = 0.0010312863960731155     (/mV2) 
  d1n = 8.397481750212752e-06     (/mV3) 
  b2n = 0.049583874211302614     (/mV) 
  c2n = -0.0003276431725574875     (/mV2) 
  d2n = 6.640950102009075e-07     (/mV3) 
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