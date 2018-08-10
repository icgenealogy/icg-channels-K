NEURON
{
  SUFFIX ch_Kdrfast 
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

  an = 0.12323818590783266     (/mV) 
  bn = -3.264328749512839     (1) 
  vhn = -29.774859653120583     (mV) 
  An = 4.411615836004768     (/ms) 
  b1n = -0.09874819868766736     (/mV) 
  c1n = 0.0010412283877650062     (/mV2) 
  d1n = -4.036267124223318e-06     (/mV3) 
  b2n = -0.03692905289062689     (/mV) 
  c2n = -0.0002915765926249723     (/mV2) 
  d2n = -2.178682660228416e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
}

STATE
{
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*n*n*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}