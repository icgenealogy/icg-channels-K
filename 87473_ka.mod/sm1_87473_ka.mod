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

  ab = -0.1999986227975416     (/mV) 
  bb = 13.999910865810211     (1) 
  vhb = -224.97933084883874     (mV) 
  Ab = 22.63565868473817     (/ms) 
  b1b = 1.622223243283904e-06     (/mV) 
  c1b = 4.7487399666839616e-05     (/mV2) 
  d1b = -7.594572945236022e-08     (/mV3) 
  b2b = 0.006123271668388627     (/mV) 
  c2b = -1.6593949213885287e-05     (/mV2) 
  d2b = 1.567784678611161e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  bInf 
  bTau 
}

STATE
{
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*b
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}