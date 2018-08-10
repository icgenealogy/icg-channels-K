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

  ap = 0.49999765367842014     (/mV) 
  bp = -17.4998863892096     (1) 
  vhp = -107.60233282821035     (mV) 
  Ap = 5.23999993281828     (/ms) 
  b1p = 1.6522230489225024e-09     (/mV) 
  c1p = 1.0211366399431137e-08     (/mV2) 
  d1p = 2.964766444824038e-11     (/mV3) 
  b2p = -4.1430213873134986e-11     (/mV) 
  c2p = 1.0240061414665531e-08     (/mV2) 
  d2p = 2.9491707002771756e-11     (/mV3) 

  aq = -0.1515142973990631     (/mV) 
  bq = 7.57575769290872     (1) 
  vhq = -54.32840868892064     (mV) 
  Aq = 49.746388388821394     (/ms) 
  b1q = 0.002479466160514819     (/mV) 
  c1q = -2.5781216419935503e-05     (/mV2) 
  d1q = 8.38797334859157e-08     (/mV3) 
  b2q = -0.10637936769283741     (/mV) 
  c2q = -0.002704520572452811     (/mV2) 
  d2q = -2.1913011815414778e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  pInf 
  pTau 
  qInf 
  qTau 
}

STATE
{
  p
  q
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*p*q
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  p' = (pInf - p) / pTau 
  q' = (qInf - q) / qTau 
}

INITIAL
{
  rates(v)
  p = pInf 
  q = qInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 

    qInf = 1/(1 + exp(-aq*v + bq)) 
    qTau = Aq / ( exp(-(b1q*(v-vhq) + c1q*(v-vhq)^2 + d1q*(v-vhq)^3)) + exp((b2q*(v-vhq) + c2q*(v-vhq)^2 + d2q*(v-vhq)^3)) ) 


  UNITSON
}