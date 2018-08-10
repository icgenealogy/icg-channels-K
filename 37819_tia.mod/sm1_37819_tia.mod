NEURON
{
  SUFFIX ia 
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

  ah = -0.16666451865443077     (/mV) 
  bh = 12.99985472647992     (1) 
  vhh = -69.99641598269255     (mV) 
  Ah = 21.16283989972838     (/ms) 
  b1h = -0.09750809763864891     (/mV) 
  c1h = 0.0016976836957911258     (/mV2) 
  d1h = -6.883440859195506e-06     (/mV3) 
  b2h = -0.012659697260573746     (/mV) 
  c2h = 0.0005070296515200169     (/mV2) 
  d2h = -2.409052794096825e-06     (/mV3) 

  am = 0.11764647781884505     (/mV) 
  bm = -7.058778104944586     (1) 
  vhm = -61.62574585343485     (mV) 
  Am = 1.066019604896712     (/ms) 
  b1m = 0.07624788134141457     (/mV) 
  c1m = 0.0005219393745997911     (/mV2) 
  d1m = -1.4490403215114496e-06     (/mV3) 
  b2m = 0.057368683791952015     (/mV) 
  c2m = -0.00044326661421633957     (/mV2) 
  d2m = 1.1364894216952166e-06     (/mV3) 
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
  g = gbar*h*m*m*m
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