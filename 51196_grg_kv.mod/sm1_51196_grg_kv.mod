NEURON
{
  SUFFIX GrG_KV 
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

  an = 0.057292485404867775     (/mV) 
  bn = -1.2021715016112837     (1) 
  vhn = -58.460630469682016     (mV) 
  An = 0.39770764258995717     (/ms) 
  b1n = 0.04075512659996585     (/mV) 
  c1n = 0.000710669365259825     (/mV2) 
  d1n = 7.517450220278197e-06     (/mV3) 
  b2n = 0.025914427124609004     (/mV) 
  c2n = -5.860894749723128e-05     (/mV2) 
  d2n = -7.73821135279106e-08     (/mV3) 
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