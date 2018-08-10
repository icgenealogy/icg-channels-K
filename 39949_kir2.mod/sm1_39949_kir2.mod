NEURON
{
  SUFFIX kir2 
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

  an = -0.09090873646880361     (/mV) 
  bn = 10.09088078135547     (1) 
  vhn = -75.00041496706008     (mV) 
  An = 0.21277022569342235     (/ms) 
  b1n = 0.0033211057169928978     (/mV) 
  c1n = -2.720518695033056e-05     (/mV2) 
  d1n = 4.5932584058202474e-07     (/mV3) 
  b2n = -0.0011996473976927248     (/mV) 
  c2n = 6.137889585459585e-05     (/mV2) 
  d2n = -1.993959706811425e-07     (/mV3) 
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
  g = gbar*n
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