NEURON
{
  SUFFIX kdrinter 
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

  an = 0.11714884007089604     (/mV) 
  bn = -2.225806068922662     (1) 
  vhn = -18.81310421842136     (mV) 
  An = 6.627848330327646     (/ms) 
  b1n = 0.02283322336891226     (/mV) 
  c1n = -9.320030829752678e-06     (/mV2) 
  d1n = -5.0183652750104244e-08     (/mV3) 
  b2n = 0.09441279456376836     (/mV) 
  c2n = -2.7869249567045168e-05     (/mV2) 
  d2n = -1.8349553227230824e-07     (/mV3) 
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