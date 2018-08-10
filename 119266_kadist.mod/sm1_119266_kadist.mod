NEURON
{
  SUFFIX kad 
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

  an = 0.07019412066677275     (/mV) 
  bn = -0.04950013705421671     (1) 
  vhn = -7.693488082999588     (mV) 
  An = 1.1491397431789536     (/ms) 
  b1n = -0.017131591098838703     (/mV) 
  c1n = -0.00019637509787629758     (/mV2) 
  d1n = 1.8833519284094368e-06     (/mV3) 
  b2n = -0.05448798231775163     (/mV) 
  c2n = -0.0002483085744757831     (/mV2) 
  d2n = 3.4552834563594747e-06     (/mV3) 

  al = -0.11222950956701772     (/mV) 
  bl = 6.284843871858036     (1) 
  vhl = -35.92917595703537     (mV) 
  Al = 7.411254539543285     (/ms) 
  b1l = 0.0760046803335803     (/mV) 
  c1l = 0.001315458376340289     (/mV2) 
  d1l = -1.0333374184705243e-05     (/mV3) 
  b2l = -0.01902154587961703     (/mV) 
  c2l = 5.4406216694296154e-05     (/mV2) 
  d2l = -3.042109101275883e-08     (/mV3) 
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
  g = gbar*n*l
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