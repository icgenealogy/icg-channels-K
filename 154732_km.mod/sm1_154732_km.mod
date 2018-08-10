NEURON
{
  SUFFIX km 
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

  an = 0.11106201810705864     (/mV) 
  bn = -3.331986447892609     (1) 
  vhn = -34.84193048439757     (mV) 
  An = 108.60256517579806     (/ms) 
  b1n = 0.0638211648920317     (/mV) 
  c1n = 0.0007865262614640948     (/mV2) 
  d1n = 4.403569818556329e-06     (/mV3) 
  b2n = 0.04516794111929679     (/mV) 
  c2n = -0.0003054668878015298     (/mV2) 
  d2n = 8.670321376820413e-07     (/mV3) 
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