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

  ah = -0.0999542345304833     (/mV) 
  bh = 8.595910838821538     (1) 
  vhh = -46.78710484507039     (mV) 
  Ah = 1.3800000009267626     (/ms) 
  b1h = 2.5681266533840525e-08     (/mV) 
  c1h = 4.302515184391077e-09     (/mV2) 
  d1h = 4.1900705207799556e-11     (/mV3) 
  b2h = 2.574225259944787e-08     (/mV) 
  c2h = 4.303710755623547e-09     (/mV2) 
  d2h = 4.1875061961925776e-11     (/mV3) 

  an = 0.10216819421912704     (/mV) 
  bn = -5.483168174286794     (1) 
  vhn = 17.62899219770324     (mV) 
  An = 0.5915712942236204     (/ms) 
  b1n = 0.00827336032010251     (/mV) 
  c1n = -0.00013686379806602853     (/mV2) 
  d1n = 7.12663336481534e-07     (/mV3) 
  b2n = -0.2819559455469909     (/mV) 
  c2n = -0.005932896005630948     (/mV2) 
  d2n = -3.233927153901667e-05     (/mV3) 
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