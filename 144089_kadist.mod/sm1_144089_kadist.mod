NEURON
{
  SUFFIX kad 
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

  an = 0.07019421003575728     (/mV) 
  bn = -0.04949861580512833     (1) 
  vhn = -23.84550280366718     (mV) 
  An = 0.7512522987630229     (/ms) 
  b1n = 0.0034461571879142367     (/mV) 
  c1n = -0.00024023702134541133     (/mV2) 
  d1n = 5.568567425587435e-07     (/mV3) 
  b2n = -0.07914885976352377     (/mV) 
  c2n = -0.00047694279977173734     (/mV2) 
  d2n = 8.14834321390785e-06     (/mV3) 

  al = -0.11222933927742808     (/mV) 
  bl = 6.284831890072787     (1) 
  vhl = -35.94810047213529     (mV) 
  Al = 7.4167471362220665     (/ms) 
  b1l = 0.018999859741188927     (/mV) 
  c1l = -5.426046965304823e-05     (/mV2) 
  d1l = 3.019828720985455e-08     (/mV3) 
  b2l = -0.0757642282215858     (/mV) 
  c2l = -0.0013050221826832807     (/mV2) 
  d2l = 1.0246480460791964e-05     (/mV3) 
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
  lInf 
  lTau 
}

STATE
{
  n
  l
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*l
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
  l' = (lInf - l) / lTau 
}

INITIAL
{
  rates(v)
  n = nInf 
  l = lInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 


  UNITSON
}