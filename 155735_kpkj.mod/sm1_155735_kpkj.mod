NEURON
{
  SUFFIX kpkj 
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

  ah = -0.026170808823382526     (/mV) 
  bh = -0.4885372262684073     (1) 
  vhh = -11.94319687305776     (mV) 
  Ah = 6.561972625156418     (/ms) 
  b1h = -0.01148761228226989     (/mV) 
  c1h = 0.00015510911864222306     (/mV2) 
  d1h = 3.0216510692033354e-06     (/mV3) 
  b2h = 0.07676972701458686     (/mV) 
  c2h = -0.001126356313621537     (/mV2) 
  d2h = 5.225890103133043e-06     (/mV3) 

  am = 0.06493479788945829     (/mV) 
  bm = -2.2727116422040528     (1) 
  vhm = -46.14553973664792     (mV) 
  Am = 8.571131965822348     (/ms) 
  b1m = -0.11063747975820046     (/mV) 
  c1m = 0.000988088487519856     (/mV2) 
  d1m = -2.9513662720119295e-06     (/mV3) 
  b2m = -0.11055676280828726     (/mV) 
  c2m = -0.002645910766323135     (/mV2) 
  d2m = -2.7424542730655438e-05     (/mV3) 
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