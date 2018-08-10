NEURON
{
  SUFFIX kap 
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

  al = -0.11223031905593263     (/mV) 
  bl = 6.284897035131566     (1) 
  vhl = -35.86675018295306     (mV) 
  Al = 7.427122105207961     (/ms) 
  b1l = 0.07583812695292327     (/mV) 
  c1l = 0.0013022391923030291     (/mV2) 
  d1l = -1.0242167609723297e-05     (/mV3) 
  b2l = -0.018992732816923673     (/mV) 
  c2l = 5.424343718875645e-05     (/mV2) 
  d2l = -3.0092734367320216e-08     (/mV3) 

  an = 0.05795829824475162     (/mV) 
  bn = 0.6520523001439968     (1) 
  vhn = -42.383933633837245     (mV) 
  An = 0.7482982883025182     (/ms) 
  b1n = 0.022415847151284394     (/mV) 
  c1n = -0.00025432239053830435     (/mV2) 
  d1n = 4.773706582953093e-08     (/mV3) 
  b2n = -0.09734976531255084     (/mV) 
  c2n = -0.0011768392300829482     (/mV2) 
  d2n = 1.2059328343421665e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  lInf 
  lTau 
  nInf 
  nTau 
}

STATE
{
  l
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*l*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  l' = (lInf - l) / lTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  l = lInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}