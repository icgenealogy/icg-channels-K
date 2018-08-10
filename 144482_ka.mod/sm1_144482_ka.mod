NEURON
{
  SUFFIX ka 
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

  ah = -0.1666496473078211     (/mV) 
  bh = 12.99872565792106     (1) 
  vhh = -69.9237582964954     (mV) 
  Ah = 46.10889808403524     (/ms) 
  b1h = 0.012752059228710314     (/mV) 
  c1h = -0.0005062078252887418     (/mV2) 
  d1h = 2.4035694090436566e-06     (/mV3) 
  b2h = 0.09658360677954222     (/mV) 
  c2h = -0.0016823971637414503     (/mV2) 
  d2h = 6.8243038235963205e-06     (/mV3) 

  am = 0.11764201726691868     (/mV) 
  bm = -7.058508137370722     (1) 
  vhm = -60.41746044077383     (mV) 
  Am = 2.3579159715306095     (/ms) 
  b1m = -0.060233985538912384     (/mV) 
  c1m = 0.0004818355750525758     (/mV2) 
  d1m = -1.2633847270118157e-06     (/mV3) 
  b2m = -0.07089901689738078     (/mV) 
  c2m = -0.00029883531859809393     (/mV2) 
  d2m = 3.922040991963368e-06     (/mV3) 
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