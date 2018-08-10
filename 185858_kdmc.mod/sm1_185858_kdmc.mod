NEURON
{
  SUFFIX kdmc 
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

  ah = -0.1249636514836142     (/mV) 
  bh = 12.497754106751758     (1) 
  vhh = -189.62286203865781     (mV) 
  Ah = 1058.6058706444753     (/ms) 
  b1h = 0.017966276440197045     (/mV) 
  c1h = -0.00012025493018309856     (/mV2) 
  d1h = 1.5008042958102844e-07     (/mV3) 
  b2h = 0.0024446053583069985     (/mV) 
  c2h = -4.5187374852276064e-05     (/mV2) 
  d2h = 1.4857524406055434e-07     (/mV3) 

  am = 0.07142846555805117     (/mV) 
  bm = -1.7857067210950008     (1) 
  vhm = -76.8469043997517     (mV) 
  Am = 2.020000000836314     (/ms) 
  b1m = 1.1239427133481175e-08     (/mV) 
  c1m = -2.151097676345664e-10     (/mV2) 
  d1m = 3.6147520747420325e-11     (/mV3) 
  b2m = 1.1186275264759563e-08     (/mV) 
  c2m = -2.1273423447890534e-10     (/mV2) 
  d2m = 3.6129611321248784e-11     (/mV3) 
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