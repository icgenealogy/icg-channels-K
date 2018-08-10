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

  al = -0.11222979612609103     (/mV) 
  bl = 6.284862656526089     (1) 
  vhl = -35.92917595703537     (mV) 
  Al = 7.411254539543285     (/ms) 
  b1l = 0.0760046803335803     (/mV) 
  c1l = 0.001315458376340289     (/mV2) 
  d1l = -1.0333374184705243e-05     (/mV3) 
  b2l = -0.01902154587961703     (/mV) 
  c2l = 5.4406216694296154e-05     (/mV2) 
  d2l = -3.042109101275883e-08     (/mV3) 

  an = 0.05795845209699409     (/mV) 
  bn = 0.652057183730017     (1) 
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