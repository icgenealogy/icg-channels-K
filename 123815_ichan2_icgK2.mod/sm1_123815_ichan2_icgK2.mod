NEURON
{
  SUFFIX ichan2 
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

  ah = -0.1458419116768305     (/mV) 
  bh = 6.935000900124219     (1) 
  vhh = -44.15471306383958     (mV) 
  Ah = 9.255808879561929     (/ms) 
  b1h = -0.12106753380406858     (/mV) 
  c1h = 0.0013362167682621363     (/mV2) 
  d1h = -4.627540004028e-06     (/mV3) 
  b2h = -0.035096254856456     (/mV) 
  c2h = 0.00019942359417236848     (/mV2) 
  d2h = 3.056680133115912e-07     (/mV3) 

  anf = 0.12312442734320939     (/mV) 
  bnf = -3.2585016968762597     (1) 
  vhnf = -29.60580001630258     (mV) 
  Anf = 6.104852829055666     (/ms) 
  b1nf = 0.03587575479571085     (/mV) 
  c1nf = 0.00025819842859158606     (/mV2) 
  d1nf = 1.8848912126251375e-06     (/mV3) 
  b2nf = 0.09858237705150741     (/mV) 
  c2nf = -0.001028712026745518     (/mV2) 
  d2nf = 3.945692898450107e-06     (/mV3) 

  am = 0.12460803098484766     (/mV) 
  bm = -3.609804207890874     (1) 
  vhm = -34.86329658157132     (mV) 
  Am = 0.2292939553968912     (/ms) 
  b1m = -0.03773332046669612     (/mV) 
  c1m = 0.00027307118759936817     (/mV2) 
  d1m = -8.110744669882462e-07     (/mV3) 
  b2m = -0.05484467321830793     (/mV) 
  c2m = -0.0008159509976556636     (/mV2) 
  d2m = -5.772887994467185e-06     (/mV3) 

  ans = 0.12327799698597372     (/mV) 
  bns = -4.745063575610305     (1) 
  vhns = -44.75691006645519     (mV) 
  Ans = 16.55066474665615     (/ms) 
  b1ns = -0.0928914700826212     (/mV) 
  c1ns = 0.0008632984510914809     (/mV2) 
  d1ns = -2.94624296969811e-06     (/mV3) 
  b2ns = -0.05070696445166038     (/mV) 
  c2ns = -0.0007145999959043583     (/mV2) 
  d2ns = -6.377119413802662e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  nfInf 
  nfTau 
  mInf 
  mTau 
  nsInf 
  nsTau 
}

STATE
{
  h
  nf
  m
  ns
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*nf*nf*nf*nf*m*m*m*ns*ns*ns*ns
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  nf' = (nfInf - nf) / nfTau 
  m' = (mInf - m) / mTau 
  ns' = (nsInf - ns) / nsTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  nf = nfInf 
  m = mInf 
  ns = nsInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nfInf = 1/(1 + exp(-anf*v + bnf)) 
    nfTau = Anf / ( exp(-(b1nf*(v-vhnf) + c1nf*(v-vhnf)^2 + d1nf*(v-vhnf)^3)) + exp((b2nf*(v-vhnf) + c2nf*(v-vhnf)^2 + d2nf*(v-vhnf)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    nsInf = 1/(1 + exp(-ans*v + bns)) 
    nsTau = Ans / ( exp(-(b1ns*(v-vhns) + c1ns*(v-vhns)^2 + d1ns*(v-vhns)^3)) + exp((b2ns*(v-vhns) + c2ns*(v-vhns)^2 + d2ns*(v-vhns)^3)) ) 


  UNITSON
}