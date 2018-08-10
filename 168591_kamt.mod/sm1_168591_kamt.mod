NEURON
{
  SUFFIX kamt 
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

  ah = -0.1666645006944623     (/mV) 
  bh = 6.9499203430458065     (1) 
  vhh = -69.99500214927993     (mV) 
  Ah = 13.331326998095667     (/ms) 
  b1h = 0.1981243577426173     (/mV) 
  c1h = -4.173336231076614e-06     (/mV2) 
  d1h = -3.383847802372381e-06     (/mV3) 
  b2h = 0.002003640297729071     (/mV) 
  c2h = 8.251673707689034e-08     (/mV2) 
  d2h = -5.734993900023243e-10     (/mV3) 

  am = 0.07142847693093747     (/mV) 
  bm = 1.25000780454829     (1) 
  vhm = -44.70739265871996     (mV) 
  Am = 6.047551266611443     (/ms) 
  b1m = 0.07410433570341611     (/mV) 
  c1m = -2.5892049496585116e-05     (/mV2) 
  d1m = -6.551562067699407e-07     (/mV3) 
  b2m = 0.025522016608467563     (/mV) 
  c2m = -6.255097622239523e-06     (/mV2) 
  d2m = 1.4851124017861437e-08     (/mV3) 
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