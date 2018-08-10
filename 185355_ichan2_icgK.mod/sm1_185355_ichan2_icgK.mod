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

  ah = -0.14583708764603898     (/mV) 
  bh = 6.934762487014717     (1) 
  vhh = -33.192235832852866     (mV) 
  Ah = 0.17520580542203681     (/ms) 
  b1h = -0.10761232418543412     (/mV) 
  c1h = 0.004281306434603029     (/mV2) 
  d1h = -2.6825318377767548e-05     (/mV3) 
  b2h = 0.014531846111444206     (/mV) 
  c2h = 0.0004935399981875977     (/mV2) 
  d2h = -3.908317628370521e-06     (/mV3) 

  anf = 0.12292552916873445     (/mV) 
  bnf = -3.2516660077017363     (1) 
  vhnf = -37.98770380747966     (mV) 
  Anf = 0.25368554382221675     (/ms) 
  b1nf = -0.08657906561993581     (/mV) 
  c1nf = 0.0011393313659135062     (/mV2) 
  d1nf = -5.190607434513382e-06     (/mV3) 
  b2nf = -0.07018285364494109     (/mV) 
  c2nf = -0.001265387149172826     (/mV2) 
  d2nf = -1.0200791737104372e-05     (/mV3) 

  am = 0.1244448170731405     (/mV) 
  bm = -3.6038692378167445     (1) 
  vhm = -33.48451615848832     (mV) 
  Am = 0.27615429298711713     (/ms) 
  b1m = 0.06685342084357203     (/mV) 
  c1m = 0.0009148183089095813     (/mV2) 
  d1m = 6.334119781807113e-06     (/mV3) 
  b2m = 0.10547295337042718     (/mV) 
  c2m = -0.0018881575745520534     (/mV2) 
  d2m = 1.25695819408574e-05     (/mV3) 

  ans = 0.12322755932563467     (/mV) 
  bns = -4.742889627264853     (1) 
  vhns = -45.325474447022984     (mV) 
  Ans = 0.5839592696636453     (/ms) 
  b1ns = 0.053320549702270195     (/mV) 
  c1ns = 0.0008200031204979145     (/mV2) 
  d1ns = 7.555821782871203e-06     (/mV3) 
  b2ns = 0.09288995599331366     (/mV) 
  c2ns = -0.0009188444084740952     (/mV2) 
  d2ns = 3.183772904054908e-06     (/mV3) 
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