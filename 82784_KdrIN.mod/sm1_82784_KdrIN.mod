NEURON
{
  SUFFIX KdrIN 
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

  an = 0.04694218319677823     (/mV) 
  bn = -1.371846079579888     (1) 
  vhn = -28.344554365513655     (mV) 
  An = 4.408864221978389     (/ms) 
  b1n = 0.010560800507169201     (/mV) 
  c1n = -2.7318883047592403e-05     (/mV2) 
  d1n = -9.061729453287502e-07     (/mV3) 
  b2n = 0.012779333976536677     (/mV) 
  c2n = 0.00012021119374041998     (/mV2) 
  d2n = -9.478271234074596e-07     (/mV3) 
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