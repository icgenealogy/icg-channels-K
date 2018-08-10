NEURON
{
  SUFFIX kl 
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

  al = -0.04999877455154249     (/mV) 
  bl = 1.5000442460607004     (1) 
  vhl = 128.94577948151672     (mV) 
  Al = 2744.705843951025     (/ms) 
  b1l = 0.03779986478335093     (/mV) 
  c1l = 0.0013448849475766614     (/mV2) 
  d1l = 7.280705533940426e-06     (/mV3) 
  b2l = -0.037798796438240315     (/mV) 
  c2l = -0.0003174750894991015     (/mV2) 
  d2l = -1.2805763390071406e-06     (/mV3) 

  an = 0.049999665109989815     (/mV) 
  bn = -2.999973965691633     (1) 
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