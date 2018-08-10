NEURON
{
  SUFFIX kv2_hh 
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

  ah = -0.026720555686588335     (/mV) 
  bh = -0.016285371650924054     (1) 
  vhh = 55.04410974516541     (mV) 
  Ah = 4984.6940620009245     (/ms) 
  b1h = 0.044978451155622254     (/mV) 
  c1h = 0.000372135503565259     (/mV2) 
  d1h = 1.024573217562434e-06     (/mV3) 
  b2h = 0.06213891050046635     (/mV) 
  c2h = -0.000869177753250626     (/mV2) 
  d2h = 5.295311447649243e-06     (/mV3) 

  an = 0.09990640125447023     (/mV) 
  bn = 1.7505445868958516     (1) 
  vhn = -30.848763299219726     (mV) 
  An = 45.156839870210526     (/ms) 
  b1n = -0.07286017795772105     (/mV) 
  c1n = 0.00078984830582902     (/mV2) 
  d1n = -2.785171016808687e-06     (/mV3) 
  b2n = -0.058730853526072484     (/mV) 
  c2n = -0.0005926280515163091     (/mV2) 
  d2n = -2.4003962035829148e-06     (/mV3) 
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
  nInf 
  nTau 
}

STATE
{
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}