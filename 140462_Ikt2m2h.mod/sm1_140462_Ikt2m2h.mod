NEURON
{
  SUFFIX Ikt2m2h 
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

  ah = -0.08128080146808837     (/mV) 
  bh = 3.2759329618800703     (1) 
  vhh = -7.161389880235354     (mV) 
  Ah = 370.44190089901144     (/ms) 
  b1h = 0.08341421318424862     (/mV) 
  c1h = 0.0018726877792682923     (/mV2) 
  d1h = 1.1700588339169356e-05     (/mV3) 
  b2h = 0.011917506711413088     (/mV) 
  c2h = 7.400829812412386e-05     (/mV2) 
  d2h = -9.59042483923325e-07     (/mV3) 

  am = 0.0758723763332122     (/mV) 
  bm = -0.7458191763501747     (1) 
  vhm = -18.44343102489529     (mV) 
  Am = 7.267380678059852     (/ms) 
  b1m = 0.0914915444940927     (/mV) 
  c1m = 0.0025930740257568757     (/mV2) 
  d1m = 1.8860176251817248e-05     (/mV3) 
  b2m = 0.0031596989032245427     (/mV) 
  c2m = 9.770309821943327e-05     (/mV2) 
  d2m = -6.559761049275117e-07     (/mV3) 
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