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

  al = -0.09087162046940205     (/mV) 
  bl = 4.907895329740452     (1) 
  vhl = -17.26000197035997     (mV) 
  Al = 463.8949738128468     (/ms) 
  b1l = -0.00894522572622225     (/mV) 
  c1l = -0.00021609449355392829     (/mV2) 
  d1l = 1.5760257427180512e-06     (/mV3) 
  b2l = 0.03184282385101712     (/mV) 
  c2l = -0.00017691950514955559     (/mV2) 
  d2l = -3.2376784061178525e-06     (/mV3) 

  an = 0.06849301324867617     (/mV) 
  bn = -0.9588967826241886     (1) 
  vhn = -60.89591529309611     (mV) 
  An = 21.944149185348973     (/ms) 
  b1n = 0.09489481256573025     (/mV) 
  c1n = 0.00286053148980949     (/mV2) 
  d1n = 3.9462520407916e-05     (/mV3) 
  b2n = 0.03369918001272974     (/mV) 
  c2n = -0.00011314778063642702     (/mV2) 
  d2n = 2.814688353885433e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  lInf 
  lTau 
  nInf 
  nTau 
}

STATE
{
  l
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*l*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  l' = (lInf - l) / lTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  l = lInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}