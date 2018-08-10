NEURON
{
  SUFFIX KdBG 
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

  axs = 0.22452251555994135     (/mV) 
  bxs = -14.144917509325074     (1) 
  vhxs = -100.00372785923228     (mV) 
  Axs = 2.020012628718731     (/ms) 
  b1xs = -0.005374252156330549     (/mV) 
  c1xs = 1.4839713448429378e-05     (/mV2) 
  d1xs = -8.990148977491942e-09     (/mV3) 
  b2xs = -0.005311621110925433     (/mV) 
  c2xs = -1.833395081284254e-05     (/mV2) 
  d2xs = 5.541750044274985e-08     (/mV3) 

  ays = -0.09285321039472165     (/mV) 
  bys = 6.786880935371288     (1) 
  vhys = -113.72949650632484     (mV) 
  Ays = 4193.437333521145     (/ms) 
  b1ys = 0.015384695028891933     (/mV) 
  c1ys = -0.0005948904618266506     (/mV2) 
  d1ys = -2.456563260384144e-06     (/mV3) 
  b2ys = -1.3338800713965615     (/mV) 
  c2ys = -0.0031171731026928677     (/mV2) 
  d2ys = 6.81729974017634e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  xsInf 
  xsTau 
  ysInf 
  ysTau 
}

STATE
{
  xs
  ys
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*xs*ys
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  xs' = (xsInf - xs) / xsTau 
  ys' = (ysInf - ys) / ysTau 
}

INITIAL
{
  rates(v)
  xs = xsInf 
  ys = ysInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    xsInf = 1/(1 + exp(-axs*v + bxs)) 
    xsTau = Axs / ( exp(-(b1xs*(v-vhxs) + c1xs*(v-vhxs)^2 + d1xs*(v-vhxs)^3)) + exp((b2xs*(v-vhxs) + c2xs*(v-vhxs)^2 + d2xs*(v-vhxs)^3)) ) 

    ysInf = 1/(1 + exp(-ays*v + bys)) 
    ysTau = Ays / ( exp(-(b1ys*(v-vhys) + c1ys*(v-vhys)^2 + d1ys*(v-vhys)^3)) + exp((b2ys*(v-vhys) + c2ys*(v-vhys)^2 + d2ys*(v-vhys)^3)) ) 


  UNITSON
}