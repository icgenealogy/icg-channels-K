NEURON
{
  SUFFIX B_A 
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

  ah = -0.09995423452390129     (/mV) 
  bh = 8.59591083825134     (1) 
  vhh = -46.78710484507039     (mV) 
  Ah = 1.3800000009267626     (/ms) 
  b1h = 2.5681266533840525e-08     (/mV) 
  c1h = 4.302515184391077e-09     (/mV2) 
  d1h = 4.1900705207799556e-11     (/mV3) 
  b2h = 2.574225259944787e-08     (/mV) 
  c2h = 4.303710755623547e-09     (/mV2) 
  d2h = 4.1875061961925776e-11     (/mV3) 

  an = 0.10216894749636292     (/mV) 
  bn = -5.483222075252121     (1) 
  vhn = 38.023125256063075     (mV) 
  An = 0.5868917420717303     (/ms) 
  b1n = 0.31027586973911564     (/mV) 
  c1n = 0.005569766048467152     (/mV2) 
  d1n = 2.5423280019765534e-05     (/mV3) 
  b2n = -0.00830070051416719     (/mV) 
  c2n = 0.00010758672743618691     (/mV2) 
  d2n = -1.962821350007798e-07     (/mV3) 
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
  g = gbar*h*n*n*n*n
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