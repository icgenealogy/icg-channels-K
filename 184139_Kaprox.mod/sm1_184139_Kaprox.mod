NEURON
{
  SUFFIX KaProx 
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

  an = 0.12051682722885787     (/mV) 
  bn = -1.7916132113059395     (1) 
  vhn = -22.870633249916253     (mV) 
  An = 2.3872044537065507     (/ms) 
  b1n = 0.08689559231370546     (/mV) 
  c1n = 8.103413483208847e-05     (/mV2) 
  d1n = -8.025422935475841e-06     (/mV3) 
  b2n = 0.03894580629757757     (/mV) 
  c2n = 0.00046601106923845867     (/mV2) 
  d2n = -5.2602005126040745e-06     (/mV3) 

  al = -0.1949968107325358     (/mV) 
  bl = 13.649798394712962     (1) 
  vhl = 175.90010707195637     (mV) 
  Al = 24.136331371573466     (/ms) 
  b1l = -0.03430680294189799     (/mV) 
  c1l = -0.0003019395805821917     (/mV2) 
  d1l = -6.159270081504664e-07     (/mV3) 
  b2l = 0.17435097310161032     (/mV) 
  c2l = -0.0006525444726839826     (/mV2) 
  d2l = 2.824066339438821e-06     (/mV3) 
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