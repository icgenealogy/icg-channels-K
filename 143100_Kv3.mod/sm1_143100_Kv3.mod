NEURON
{
  SUFFIX Kv3 
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

  ah = -0.01364135542893582     (/mV) 
  bh = -1.2493717253875696     (1) 
  vhh = -0.0019586982287523322     (mV) 
  Ah = 33.75339361867031     (/ms) 
  b1h = -0.039958744630089985     (/mV) 
  c1h = -8.656362196573234e-05     (/mV2) 
  d1h = 4.6435931928430915e-06     (/mV3) 
  b2h = -0.03996374053792239     (/mV) 
  c2h = 8.651039149661527e-05     (/mV2) 
  d2h = 4.644090428145206e-06     (/mV3) 

  am = 0.12820467523344262     (/mV) 
  bm = -3.3333128789651156     (1) 
  vhm = -24.897754987081875     (mV) 
  Am = 14.02296765213098     (/ms) 
  b1m = -0.08875792549011133     (/mV) 
  c1m = 0.00011338189710213709     (/mV2) 
  d1m = 2.768831921774766e-06     (/mV3) 
  b2m = -0.0686925759208971     (/mV) 
  c2m = 0.00025642481725871914     (/mV2) 
  d2m = 4.758986236759531e-06     (/mV3) 
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