NEURON
{
  SUFFIX K_Pst 
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

  ah = -0.09086598695861259     (/mV) 
  bh = 5.816404858189466     (1) 
  vhh = -28.555489121158317     (mV) 
  Ah = 561.5480819646447     (/ms) 
  b1h = -0.006459228276704687     (/mV) 
  c1h = -0.00021572512470089411     (/mV2) 
  d1h = 1.3628403504112846e-06     (/mV3) 
  b2h = 0.03380040414939692     (/mV) 
  c2h = -0.00021801013032355434     (/mV2) 
  d2h = -2.170517085630496e-06     (/mV3) 

  am = 0.08333323044071445     (/mV) 
  bm = -0.9166579541377337     (1) 
  vhm = -71.59701085183103     (mV) 
  Am = 24.9931521430333     (/ms) 
  b1m = -0.030859758709265238     (/mV) 
  c1m = 7.004755111204131e-05     (/mV2) 
  d1m = -9.387931646120578e-08     (/mV3) 
  b2m = -0.09657596422622643     (/mV) 
  c2m = -0.0038130501456876548     (/mV2) 
  d2m = -7.297363313864291e-05     (/mV3) 
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