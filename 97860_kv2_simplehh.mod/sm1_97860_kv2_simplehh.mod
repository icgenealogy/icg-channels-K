NEURON
{
  SUFFIX kv2_simplehh 
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

  an = 0.12453619402922904     (/mV) 
  bn = 0.007936135706137663     (1) 
  vhn = -30.812981493844976     (mV) 
  An = 45.14728526360964     (/ms) 
  b1n = -0.07294251544600112     (/mV) 
  c1n = 0.0007917949018482795     (/mV2) 
  d1n = -2.7955959688055603e-06     (/mV3) 
  b2n = -0.05870664013472017     (/mV) 
  c2n = -0.0005928737672171843     (/mV2) 
  d2n = -2.4061267860346243e-06     (/mV3) 
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