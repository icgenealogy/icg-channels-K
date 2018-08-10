NEURON
{
  SUFFIX Kv1 
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

  an = 0.10999804437295606     (/mV) 
  bn = -4.949897333029179     (1) 
  vhn = -44.805118771700855     (mV) 
  An = 1.5114032946786897     (/ms) 
  b1n = 0.07744030901903773     (/mV) 
  c1n = -0.00014938416647485708     (/mV2) 
  d1n = -2.282244600927644e-06     (/mV3) 
  b2n = 0.028711781390817287     (/mV) 
  c2n = 2.310795460128703e-05     (/mV2) 
  d2n = -2.0109139789987605e-07     (/mV3) 
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