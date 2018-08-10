NEURON
{
  SUFFIX Kv4 
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

  ah = -0.10403260498837612     (/mV) 
  bh = 7.186474495727684     (1) 
  vhh = -25.529546100796523     (mV) 
  Ah = 10.904229516228918     (/ms) 
  b1h = -0.008328339166849524     (/mV) 
  c1h = -1.4489315149694174e-05     (/mV2) 
  d1h = 3.617538532520666e-07     (/mV3) 
  b2h = 0.00832767823449958     (/mV) 
  c2h = -0.0004900199743788192     (/mV2) 
  d2h = 2.933220947595659e-06     (/mV3) 

  an = 0.05771312860095816     (/mV) 
  bn = -3.2896469930821826     (1) 
  vhn = -75.15894377585313     (mV) 
  An = 1.3136951144835072     (/ms) 
  b1n = -0.014074054695980063     (/mV) 
  c1n = -0.0001515015469451284     (/mV2) 
  d1n = 5.176409645450686e-07     (/mV3) 
  b2n = -0.0369079860981973     (/mV) 
  c2n = -0.00019061814844850667     (/mV2) 
  d2n = -3.709276169767865e-06     (/mV3) 
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
  nInf 
  nTau 
}

STATE
{
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*h*n*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}