NEURON
{
  SUFFIX kargc 
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

  ah = -0.1574729015581507     (/mV) 
  bh = 9.763379777219916     (1) 
  vhh = -65.0566532659131     (mV) 
  Ah = 50.00216765635073     (/ms) 
  b1h = -0.00022551977736713286     (/mV) 
  c1h = 1.3363881445955168e-05     (/mV2) 
  d1h = -6.951268601244225e-08     (/mV3) 
  b2h = -0.00021284159771237633     (/mV) 
  c2h = 1.2609076962974631e-05     (/mV2) 
  d2h = -6.5577323434682e-08     (/mV3) 

  am = 0.09958418580401021     (/mV) 
  bm = -2.1733767450222508     (1) 
  vhm = -26.88233502435725     (mV) 
  Am = 10.32114402589684     (/ms) 
  b1m = -0.06419631765899553     (/mV) 
  c1m = 0.0005309861396135052     (/mV2) 
  d1m = -1.764032721730337e-06     (/mV3) 
  b2m = -0.046865555817896885     (/mV) 
  c2m = -0.00029692047263868925     (/mV2) 
  d2m = -2.0746189996634623e-06     (/mV3) 
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
  g = gbar*h*m*m*m
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