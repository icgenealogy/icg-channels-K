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

  al = -0.1219508703723467     (/mV) 
  bl = 7.073152224683065     (1) 
  vhl = -34.91889608729459     (mV) 
  Al = 9.684725590771196     (/ms) 
  b1l = 0.006273492403452145     (/mV) 
  c1l = 8.074501741984894e-05     (/mV2) 
  d1l = -4.0829108780805385e-07     (/mV3) 
  b2l = -0.006274350906568199     (/mV) 
  c2l = -0.000664331144627411     (/mV2) 
  d2l = -8.972725957679777e-06     (/mV3) 

  an = 0.028571371669477725     (/mV) 
  bn = -0.6085683812053577     (1) 
  vhn = -107.85265658954941     (mV) 
  An = 0.37127382325521063     (/ms) 
  b1n = -0.012028964404348676     (/mV) 
  c1n = 7.60814580794719e-05     (/mV2) 
  d1n = -1.5697322366129048e-07     (/mV3) 
  b2n = -0.05815958096172296     (/mV) 
  c2n = 0.001614039351232231     (/mV2) 
  d2n = -1.836136353254895e-05     (/mV3) 
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