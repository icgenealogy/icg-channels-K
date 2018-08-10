NEURON
{
  SUFFIX Golgi_KV 
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

  an = 0.05539631723438578     (/mV) 
  bn = -1.204907965174347     (1) 
  vhn = -8.792131836452993     (mV) 
  An = 0.25206606354369165     (/ms) 
  b1n = -0.04117815386124939     (/mV) 
  c1n = 0.000541054807855945     (/mV2) 
  d1n = -1.628224359326423e-06     (/mV3) 
  b2n = -0.00010880938307231723     (/mV) 
  c2n = 7.805047095948166e-05     (/mV2) 
  d2n = 2.6621684433758616e-09     (/mV3) 
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