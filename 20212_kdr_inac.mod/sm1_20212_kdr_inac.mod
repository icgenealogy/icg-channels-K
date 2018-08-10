NEURON
{
  SUFFIX kdr_inac 
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

  an = 0.3328855270113575     (/mV) 
  bn = -13.313464930504468     (1) 
  vhn = -27.432150885145276     (mV) 
  An = 5.994800486746047     (/ms) 
  b1n = 0.01572758531508908     (/mV) 
  c1n = 0.00013063016948274968     (/mV2) 
  d1n = 4.0670969238506527e-07     (/mV3) 
  b2n = 0.015547357116596525     (/mV) 
  c2n = -0.00011933339140141391     (/mV2) 
  d2n = 3.1198520869970005e-07     (/mV3) 

  al = -0.017572572789004454     (/mV) 
  bl = 0.2546879357552057     (1) 
  vhl = -172.02866837435135     (mV) 
  Al = 99.49879342239477     (/ms) 
  b1l = 0.00011157446602529398     (/mV) 
  c1l = -6.904996604027541e-07     (/mV2) 
  d1l = 1.2990288127040785e-09     (/mV3) 
  b2l = -37.399703652299316     (/mV) 
  c2l = 0.6475355277069418     (/mV2) 
  d2l = -0.0028045425188191436     (/mV3) 
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
  lInf 
  lTau 
}

STATE
{
  n
  l
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*n*l
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
  l' = (lInf - l) / lTau 
}

INITIAL
{
  rates(v)
  n = nInf 
  l = lInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 


  UNITSON
}