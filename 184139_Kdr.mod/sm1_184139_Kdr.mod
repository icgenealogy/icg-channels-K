NEURON
{
  SUFFIX Kdr 
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

  an = 0.11699560140263643     (/mV) 
  bn = -0.4445383959876926     (1) 
  vhn = -4.3125613913227     (mV) 
  An = 50.44283802363351     (/ms) 
  b1n = -0.054284520856064     (/mV) 
  c1n = -0.0002827881647314877     (/mV2) 
  d1n = 5.463430528471535e-06     (/mV3) 
  b2n = -0.05772295188280416     (/mV) 
  c2n = 0.00022881815059978296     (/mV2) 
  d2n = 5.464210807085608e-06     (/mV3) 
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