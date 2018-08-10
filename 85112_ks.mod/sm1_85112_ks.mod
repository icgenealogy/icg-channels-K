NEURON
{
  SUFFIX ks 
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

  an = 0.053990554840282805     (/mV) 
  bn = -1.9137022977604248     (1) 
  vhn = -33.22382085859634     (mV) 
  An = 106.46261500415024     (/ms) 
  b1n = 0.02656258266980287     (/mV) 
  c1n = 9.792778481332372e-05     (/mV2) 
  d1n = 1.4598916904562265e-07     (/mV3) 
  b2n = 0.02994532908232717     (/mV) 
  c2n = -0.00010721793386873     (/mV2) 
  d2n = 1.6454223784146757e-07     (/mV3) 
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