NEURON
{
  SUFFIX Kdr 
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

  ah = -0.24919169435159305     (/mV) 
  bh = 6.232633826159075     (1) 
  vhh = -22.414295664683582     (mV) 
  Ah = 1204.2756680727985     (/ms) 
  b1h = -2.0263334633585455     (/mV) 
  c1h = 0.03521706568335249     (/mV2) 
  d1h = -0.00015605971955310328     (/mV3) 
  b2h = -0.0003046972998710284     (/mV) 
  c2h = -7.278860947365531e-06     (/mV2) 
  d2h = -5.164487796362649e-08     (/mV3) 

  am = 0.08486617714708304     (/mV) 
  bm = -0.985388399433121     (1) 
  vhm = -29.072872863127373     (mV) 
  Am = 8.876425201079556     (/ms) 
  b1m = -0.05658745196531547     (/mV) 
  c1m = 0.00040146671187544185     (/mV2) 
  d1m = -1.1822261575333273e-06     (/mV3) 
  b2m = -0.027844924582967737     (/mV) 
  c2m = 7.270788058555954e-05     (/mV2) 
  d2m = 3.232816359618065e-07     (/mV3) 
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