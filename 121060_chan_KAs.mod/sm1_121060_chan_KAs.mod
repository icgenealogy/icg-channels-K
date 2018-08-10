NEURON
{
  SUFFIX KAs 
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

  ah = -0.04597282696724026     (/mV) 
  bh = 1.5390021584344187     (1) 
  vhh = 21.875960596057713     (mV) 
  Ah = 699.8794205932867     (/ms) 
  b1h = -0.010261687237601937     (/mV) 
  c1h = 0.00016555014640240025     (/mV2) 
  d1h = 1.94196305961689e-06     (/mV3) 
  b2h = 0.006990056339712932     (/mV) 
  c2h = 2.866160191067116e-05     (/mV2) 
  d2h = 8.620050149230417e-09     (/mV3) 

  am = 0.06232608820027832     (/mV) 
  bm = -1.6821620776657362     (1) 
  vhm = -16.553357683028292     (mV) 
  Am = 14.73150405354869     (/ms) 
  b1m = -0.07599767681364306     (/mV) 
  c1m = -0.00026206383533688987     (/mV2) 
  d1m = 8.477057196746326e-06     (/mV3) 
  b2m = 0.0024070838703613836     (/mV) 
  c2m = 0.0004956599681578823     (/mV2) 
  d2m = -1.983050923358664e-06     (/mV3) 
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
  g = gbar*h*m*m
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