NEURON
{
  SUFFIX k2 
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

  ah = -0.09432443927279743     (/mV) 
  bh = 5.47088840091759     (1) 
  vhh = -229.75118679018038     (mV) 
  Ah = 79.0758397401707     (/ms) 
  b1h = -0.0026289862608230557     (/mV) 
  c1h = 8.555926790379328e-06     (/mV2) 
  d1h = -9.274384094534989e-09     (/mV3) 
  b2h = -0.035537342576295206     (/mV) 
  c2h = 0.00017681449551422995     (/mV2) 
  d2h = -5.6031478869642e-07     (/mV3) 

  am = 0.05734758432860712     (/mV) 
  bm = -0.5658459676231674     (1) 
  vhm = -39.451096628554986     (mV) 
  Am = 77.1985163659754     (/ms) 
  b1m = 0.04464299515275333     (/mV) 
  c1m = -3.747201384470792e-05     (/mV2) 
  d1m = -2.470751467660003e-06     (/mV3) 
  b2m = 0.04711477510480548     (/mV) 
  c2m = -0.00023173143378723556     (/mV2) 
  d2m = -2.1773234213640444e-08     (/mV3) 
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
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}