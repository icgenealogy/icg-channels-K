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

  ah = -0.20407040062164528     (/mV) 
  bh = 11.611662808175485     (1) 
  vhh = 38.67601038408223     (mV) 
  Ah = 21.74831339682113     (/ms) 
  b1h = -0.004986185491279878     (/mV) 
  c1h = -8.832447406869742e-06     (/mV2) 
  d1h = 4.349666359719195e-09     (/mV3) 
  b2h = 0.004986545102257126     (/mV) 
  c2h = -0.0001587488653021456     (/mV2) 
  d2h = 6.969549378350207e-07     (/mV3) 

  am = 0.11494239424502999     (/mV) 
  bm = -3.1264239014385358     (1) 
  vhm = -15.229878857732928     (mV) 
  Am = 7.336796332804628     (/ms) 
  b1m = -0.05632372438724387     (/mV) 
  c1m = 0.0007070563693287318     (/mV2) 
  d1m = -2.4742245710515515e-06     (/mV3) 
  b2m = 0.009594764626743099     (/mV) 
  c2m = 3.6601028422451435e-05     (/mV2) 
  d2m = -1.730009308336635e-07     (/mV3) 
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