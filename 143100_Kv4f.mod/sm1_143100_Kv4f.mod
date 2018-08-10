NEURON
{
  SUFFIX Kv4f 
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

  ah = -0.09999920867740662     (/mV) 
  bh = 8.299950625216105     (1) 
  vhh = -86.05475228747825     (mV) 
  Ah = 26.775156097274163     (/ms) 
  b1h = -0.05060701551788502     (/mV) 
  c1h = 0.0006096525034341445     (/mV2) 
  d1h = -1.999400378571787e-06     (/mV3) 
  b2h = -0.053574468902329184     (/mV) 
  c2h = 0.0008485896468445205     (/mV2) 
  d2h = -3.0615980229329186e-06     (/mV3) 

  am = 0.07999985873826491     (/mV) 
  bm = -3.919991211902576     (1) 
  vhm = -42.388456324920995     (mV) 
  Am = 7.089553962189055     (/ms) 
  b1m = 0.024955716469939323     (/mV) 
  c1m = -0.00010864425562009888     (/mV2) 
  d1m = -9.376356185428481e-07     (/mV3) 
  b2m = 0.03918334736088723     (/mV) 
  c2m = -8.735828104447603e-05     (/mV2) 
  d2m = -2.7063654001446106e-07     (/mV3) 
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