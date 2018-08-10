NEURON
{
  SUFFIX kd 
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

  an = 0.07449754170767386     (/mV) 
  bn = -3.4478473534171488     (1) 
  vhn = -35.11246915654801     (mV) 
  An = 2.9869057531958014     (/ms) 
  b1n = -0.06297523662308503     (/mV) 
  c1n = 0.00037170344698363774     (/mV2) 
  d1n = -7.773292572856213e-07     (/mV3) 
  b2n = -0.0036312790420648224     (/mV) 
  c2n = 9.246072818927044e-05     (/mV2) 
  d2n = 5.608269890142536e-07     (/mV3) 
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