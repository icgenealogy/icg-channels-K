NEURON
{
  SUFFIX kv 
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

  an = 0.1110597453295364     (/mV) 
  bn = 0.4750417996818126     (1) 
  vhn = 8.91202963717225     (mV) 
  An = 6.040255538821663     (/ms) 
  b1n = 0.024070510671400373     (/mV) 
  c1n = 0.00011676109876263647     (/mV2) 
  d1n = 2.6216613512602583e-07     (/mV3) 
  b2n = 0.08359350106349589     (/mV) 
  c2n = -0.0007973582491385173     (/mV2) 
  d2n = 3.1265397066897337e-06     (/mV3) 
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