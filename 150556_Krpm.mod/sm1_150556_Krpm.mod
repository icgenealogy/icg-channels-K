NEURON
{
  SUFFIX Krpm 
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

  ah = -0.05253057018626714     (/mV) 
  bh = 2.8887672577235644     (1) 
  vhh = -39.989514374275444     (mV) 
  Ah = 2456.2789545361675     (/ms) 
  b1h = 0.12489544897020063     (/mV) 
  c1h = 0.0021456495576619578     (/mV2) 
  d1h = -2.133272482552228e-05     (/mV3) 
  b2h = 0.0015647570980662619     (/mV) 
  c2h = 0.00016521725248172308     (/mV2) 
  d2h = -1.1216494816978129e-06     (/mV3) 

  am = 0.08264453638063912     (/mV) 
  bm = -1.1074303099999068     (1) 
  vhm = -54.24000663356989     (mV) 
  Am = 52.17421359414252     (/ms) 
  b1m = -0.037239047627547775     (/mV) 
  c1m = -6.52540932055061e-06     (/mV2) 
  d1m = 3.0751386803182726e-08     (/mV3) 
  b2m = -0.038188125825793146     (/mV) 
  c2m = -5.350238414512061e-06     (/mV2) 
  d2m = -1.3890298915263635e-08     (/mV3) 
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