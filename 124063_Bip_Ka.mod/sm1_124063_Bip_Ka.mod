NEURON
{
  SUFFIX IA 
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

  ah = -0.1304954744072597     (/mV) 
  bh = 9.57580774799974     (1) 
  vhh = -65.3370866391208     (mV) 
  Ah = 0.08036747568561048     (/ms) 
  b1h = -0.05803141587011992     (/mV) 
  c1h = 0.0008209819897990842     (/mV2) 
  d1h = -3.075336533755456e-06     (/mV3) 
  b2h = -0.0307967948954124     (/mV) 
  c2h = 0.0006520651919091254     (/mV2) 
  d2h = -2.7604414347523178e-06     (/mV3) 

  am = 0.13334175692399655     (/mV) 
  bm = -3.3840405312561073     (1) 
  vhm = -59.90825290677018     (mV) 
  Am = 0.05043876590167478     (/ms) 
  b1m = 4.608159380533623     (/mV) 
  c1m = -0.12351308787834296     (/mV2) 
  d1m = 0.0006540807575106251     (/mV3) 
  b2m = 0.06393164940383221     (/mV) 
  c2m = -0.0012138484917929135     (/mV2) 
  d2m = 6.012226875410461e-06     (/mV3) 
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