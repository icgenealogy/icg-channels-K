NEURON
{
  SUFFIX M99Ka 
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

  ah = -0.10999977046947292     (/mV) 
  bh = 6.159985406035528     (1) 
  vhh = -35.86675018295306     (mV) 
  Ah = 7.427122105207961     (/ms) 
  b1h = 0.07583812695292327     (/mV) 
  c1h = 0.0013022391923030291     (/mV2) 
  d1h = -1.0242167609723297e-05     (/mV3) 
  b2h = -0.018992732816923673     (/mV) 
  c2h = 5.424343718875645e-05     (/mV2) 
  d2h = -3.0092734367320216e-08     (/mV3) 

  am = 0.05878676451936042     (/mV) 
  bm = 0.6604100659327464     (1) 
  vhm = 37.71270604669929     (mV) 
  Am = 2.8542141073036382     (/ms) 
  b1m = -0.04187820126813687     (/mV) 
  c1m = 0.00017005874493832457     (/mV2) 
  d1m = -9.533523189113493e-07     (/mV3) 
  b2m = 0.00040081218825687375     (/mV) 
  c2m = 0.00021546343907177158     (/mV2) 
  d2m = 5.878536713152347e-07     (/mV3) 
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