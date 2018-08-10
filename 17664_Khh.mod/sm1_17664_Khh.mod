NEURON
{
  SUFFIX Khh 
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

  an = 0.056402884161287234     (/mV) 
  bn = -2.8788537443879085     (1) 
  vhn = -78.91575870182315     (mV) 
  An = 11.577660085062252     (/ms) 
  b1n = 0.03895364210015931     (/mV) 
  c1n = 0.0005702532325476402     (/mV2) 
  d1n = 6.255022235746779e-06     (/mV3) 
  b2n = 0.03645989120323775     (/mV) 
  c2n = -0.00017758563796952263     (/mV2) 
  d2n = 3.602704385193779e-07     (/mV3) 
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