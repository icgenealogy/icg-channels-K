NEURON
{
  SUFFIX HHk 
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

  an = 0.03161656834540454     (/mV) 
  bn = 0.14878911670629383     (1) 
  vhn = -115.09885492530688     (mV) 
  An = 1.4811362488899023     (/ms) 
  b1n = 0.00531633304269434     (/mV) 
  c1n = -5.439266141341474e-05     (/mV2) 
  d1n = 9.683478687957701e-08     (/mV3) 
  b2n = -0.018513423905139857     (/mV) 
  c2n = -0.00011340979144842395     (/mV2) 
  d2n = -1.0642598368663357e-06     (/mV3) 
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