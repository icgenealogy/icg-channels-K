NEURON
{
  SUFFIX Kv3 
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

  an = 0.07547104136312745     (/mV) 
  bn = -1.2075243754307154     (1) 
  vhn = -18.35122276327444     (mV) 
  An = 1.03035008024603     (/ms) 
  b1n = 0.039414619822909724     (/mV) 
  c1n = -4.993515691566233e-06     (/mV2) 
  d1n = -4.479071630505443e-07     (/mV3) 
  b2n = 0.03255733579973857     (/mV) 
  c2n = 9.312452749686133e-05     (/mV2) 
  d2n = -6.4770096984185e-07     (/mV3) 
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
  g = gbar*n*n
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