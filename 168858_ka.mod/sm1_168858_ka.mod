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

  ah = -0.1666602924737128     (/mV) 
  bh = 12.999563138552961     (1) 
  vhh = -69.9237582964954     (mV) 
  Ah = 46.10889808403524     (/ms) 
  b1h = 0.012752059228710314     (/mV) 
  c1h = -0.0005062078252887418     (/mV2) 
  d1h = 2.4035694090436566e-06     (/mV3) 
  b2h = 0.09658360677954222     (/mV) 
  c2h = -0.0016823971637414503     (/mV2) 
  d2h = 6.8243038235963205e-06     (/mV3) 

  am = 0.11764647781884005     (/mV) 
  bm = -7.058778104944288     (1) 
  vhm = -60.67643255034253     (mV) 
  Am = 2.3549870641179678     (/ms) 
  b1m = -0.059671819446983625     (/mV) 
  c1m = 0.0004732710526316983     (/mV2) 
  d1m = -1.2472265235636963e-06     (/mV3) 
  b2m = -0.07195941647852114     (/mV) 
  c2m = -0.0003250803030413409     (/mV2) 
  d2m = 3.876960289219264e-06     (/mV3) 
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