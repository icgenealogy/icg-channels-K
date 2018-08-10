NEURON
{
  SUFFIX KM 
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

  au = 0.09666952546766434     (/mV) 
  bu = -5.094481471268137     (1) 
  vhu = -53.071864972773774     (mV) 
  Au = 6.56133854197147     (/ms) 
  b1u = -0.041810270025777196     (/mV) 
  c1u = -3.8461703454682366e-05     (/mV2) 
  d1u = 3.007397908859068e-07     (/mV3) 
  b2u = -0.05309145868746031     (/mV) 
  c2u = 3.3550509044010776e-05     (/mV2) 
  d2u = 6.299605513747448e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  uInf 
  uTau 
}

STATE
{
  u
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*u*u
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  u' = (uInf - u) / uTau 
}

INITIAL
{
  rates(v)
  u = uInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    uInf = 1/(1 + exp(-au*v + bu)) 
    uTau = Au / ( exp(-(b1u*(v-vhu) + c1u*(v-vhu)^2 + d1u*(v-vhu)^3)) + exp((b2u*(v-vhu) + c2u*(v-vhu)^2 + d2u*(v-vhu)^3)) ) 


  UNITSON
}