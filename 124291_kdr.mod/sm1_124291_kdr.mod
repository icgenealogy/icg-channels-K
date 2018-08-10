NEURON
{
  SUFFIX kdr 
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

  an = 0.011280523115846059     (/mV) 
  bn = 0.26474223065967323     (1) 
  vhn = -99.99706835580719     (mV) 
  An = 51.40382979076925     (/ms) 
  b1n = -8.827966870247451     (/mV) 
  c1n = 0.0013403127304427744     (/mV2) 
  d1n = -4.499339074929844e-06     (/mV3) 
  b2n = -8.316190516007648     (/mV) 
  c2n = -0.002431621758876434     (/mV2) 
  d2n = 4.480362086927461e-06     (/mV3) 
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