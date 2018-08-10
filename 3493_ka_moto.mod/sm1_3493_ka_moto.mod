NEURON
{
  SUFFIX kamot 
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

  ah = -0.07575749860783945     (/mV) 
  bh = 6.500004772645109     (1) 
  vhh = -69.49195415347768     (mV) 
  Ah = 26.26660247879289     (/ms) 
  b1h = 0.02002493071635501     (/mV) 
  c1h = -0.0007784689192212638     (/mV2) 
  d1h = 3.3695605706452996e-06     (/mV3) 
  b2h = 0.044009966514418124     (/mV) 
  c2h = -0.0002698713188866922     (/mV2) 
  d2h = 1.1643983663671653e-06     (/mV3) 

  am = 0.10988994429758164     (/mV) 
  bm = -3.9780095693056374     (1) 
  vhm = -36.75523683693634     (mV) 
  Am = 1.825804192783763     (/ms) 
  b1m = -0.03255610981575679     (/mV) 
  c1m = 0.00018899257374246404     (/mV2) 
  d1m = -8.309184332308732e-07     (/mV3) 
  b2m = 0.004532099539152203     (/mV) 
  c2m = 0.0002979779294094877     (/mV2) 
  d2m = -1.3077983603214452e-06     (/mV3) 
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