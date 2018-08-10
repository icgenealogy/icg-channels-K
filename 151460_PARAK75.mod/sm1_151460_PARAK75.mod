NEURON
{
  SUFFIX parak75 
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

  an = 0.1151403896064036     (/mV) 
  bn = -7.929294715630831     (1) 
  vhn = -70.6414295498638     (mV) 
  An = 1.0918962058431534     (/ms) 
  b1n = -0.028435957113009325     (/mV) 
  c1n = 0.00012824927193131426     (/mV2) 
  d1n = -2.421161903752867e-07     (/mV3) 
  b2n = -0.019831522676185775     (/mV) 
  c2n = 0.000388648342387962     (/mV2) 
  d2n = -2.924258361887593e-06     (/mV3) 
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