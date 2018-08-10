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

  an = 0.16424049848604666     (/mV) 
  bn = -7.958220984521892     (1) 
  vhn = -49.028956272346626     (mV) 
  An = 1.7913313898643168     (/ms) 
  b1n = -0.02681405218194495     (/mV) 
  c1n = -3.9949663307814084e-05     (/mV2) 
  d1n = 2.254948087176517e-07     (/mV3) 
  b2n = -0.139523193026829     (/mV) 
  c2n = -0.00019778899978897997     (/mV2) 
  d2n = 1.7273698362746005e-06     (/mV3) 
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