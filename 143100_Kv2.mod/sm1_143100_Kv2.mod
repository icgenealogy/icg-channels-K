NEURON
{
  SUFFIX Kv2 
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

  ah = -0.0401480905464897     (/mV) 
  bh = 0.2562792843528164     (1) 
  vhh = -12.636421907180207     (mV) 
  Ah = 6800.01181642426     (/ms) 
  b1h = -0.0023396248966882183     (/mV) 
  c1h = 2.5857130951857444e-06     (/mV2) 
  d1h = -7.142490187402879e-10     (/mV3) 
  b2h = -0.0023397337330743405     (/mV) 
  c2h = -2.888531137007179e-06     (/mV2) 
  d2h = -1.399178969545917e-09     (/mV3) 

  am = 0.10988988685729198     (/mV) 
  bm = -3.6483387081192866     (1) 
  vhm = -28.40262761848712     (mV) 
  Am = 2.835633822863769     (/ms) 
  b1m = 0.032188944701286144     (/mV) 
  c1m = -0.00016441504879497235     (/mV2) 
  d1m = -1.7621013174700347e-06     (/mV3) 
  b2m = 0.08830472958225213     (/mV) 
  c2m = -0.0007239754529911368     (/mV2) 
  d2m = 1.6296066961950496e-06     (/mV3) 
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
  g = gbar*h*m*m*m*m
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