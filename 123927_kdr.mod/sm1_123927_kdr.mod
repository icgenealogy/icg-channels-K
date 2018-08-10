NEURON
{
  SUFFIX kdrG 
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

  an = 0.09999880160709167     (/mV) 
  bn = -3.499943623411191     (1) 
  vhn = -19.397929885826834     (mV) 
  An = 0.5848828911943956     (/ms) 
  b1n = -0.056085700547210925     (/mV) 
  c1n = 0.0006377326412429185     (/mV2) 
  d1n = -2.4000429647586433e-06     (/mV3) 
  b2n = -0.0036872113854201287     (/mV) 
  c2n = -0.000110072720806549     (/mV2) 
  d2n = -7.523089241081879e-07     (/mV3) 
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