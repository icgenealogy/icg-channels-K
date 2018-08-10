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

  an = 0.07019374448134359     (/mV) 
  bn = -0.04950921274894109     (1) 
  vhn = -6.545703804143589     (mV) 
  An = 1.1714800923661328     (/ms) 
  b1n = 0.053091175027905696     (/mV) 
  c1n = 0.00020603077118038567     (/mV2) 
  d1n = -3.557954293820871e-06     (/mV3) 
  b2n = 0.018470359648521596     (/mV) 
  c2n = 0.0001926088470437     (/mV2) 
  d2n = -1.991430972059963e-06     (/mV3) 

  al = -0.11223025470131802     (/mV) 
  bl = 6.284893020392546     (1) 
  vhl = -35.86675018295306     (mV) 
  Al = 7.427122105207961     (/ms) 
  b1l = 0.07583812695292327     (/mV) 
  c1l = 0.0013022391923030291     (/mV2) 
  d1l = -1.0242167609723297e-05     (/mV3) 
  b2l = -0.018992732816923673     (/mV) 
  c2l = 5.424343718875645e-05     (/mV2) 
  d2l = -3.0092734367320216e-08     (/mV3) 
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