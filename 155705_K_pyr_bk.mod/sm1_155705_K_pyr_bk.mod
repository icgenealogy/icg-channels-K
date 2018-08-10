NEURON
{
  SUFFIX KPyr 
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

  an = 0.09067052750490233     (/mV) 
  bn = -3.6990799761426723     (1) 
  vhn = -15.098073029714818     (mV) 
  An = 1.451263314480069     (/ms) 
  b1n = -0.013334405862201209     (/mV) 
  c1n = -0.00016652720847193833     (/mV2) 
  d1n = 1.3303047927599065e-06     (/mV3) 
  b2n = 0.0334129105495324     (/mV) 
  c2n = -0.0007777481981231062     (/mV2) 
  d2n = 4.612593365401464e-06     (/mV3) 
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