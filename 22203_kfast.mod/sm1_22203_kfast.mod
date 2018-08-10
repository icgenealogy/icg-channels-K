NEURON
{
  SUFFIX kfast 
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

  aa = 0.03832253419331342     (/mV) 
  ba = -1.6737997964311555     (1) 
  vha = -4.467073187529895     (mV) 
  Aa = 1.2716417633604742     (/ms) 
  b1a = -0.0030952354519848857     (/mV) 
  c1a = -4.3045691462039304e-05     (/mV2) 
  d1a = 2.277950096912976e-08     (/mV3) 
  b2a = 0.02506742434338708     (/mV) 
  c2a = -0.0002117848849763933     (/mV2) 
  d2a = 2.667426504060322e-07     (/mV3) 

  ab = -0.1393663571833226     (/mV) 
  bb = 10.280963926409166     (1) 
  vhb = -78.20200212541535     (mV) 
  Ab = 123.49821240668967     (/ms) 
  b1b = -0.09162849114067927     (/mV) 
  c1b = 0.0009032209101579718     (/mV2) 
  d1b = -2.7158144653118493e-06     (/mV3) 
  b2b = -0.12842005243211652     (/mV) 
  c2b = -0.004906532860299506     (/mV2) 
  d2b = -0.0001152996566071113     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  aInf 
  aTau 
  bInf 
  bTau 
}

STATE
{
  a
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*a*a*a*a*b
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  a' = (aInf - a) / aTau 
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  a = aInf 
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    aInf = 1/(1 + exp(-aa*v + ba)) 
    aTau = Aa / ( exp(-(b1a*(v-vha) + c1a*(v-vha)^2 + d1a*(v-vha)^3)) + exp((b2a*(v-vha) + c2a*(v-vha)^2 + d2a*(v-vha)^3)) ) 

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}