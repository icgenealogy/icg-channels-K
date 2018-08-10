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

  al = -0.11222973351354619     (/mV) 
  bl = 6.2848578408636495     (1) 
  vhl = -35.92112876558254     (mV) 
  Al = 7.416245205560595     (/ms) 
  b1l = 0.07600259011307714     (/mV) 
  c1l = 0.0013100157485364007     (/mV2) 
  d1l = -1.0295401214981123e-05     (/mV3) 
  b2l = -0.01900783998622462     (/mV) 
  c2l = 5.4305735352585535e-05     (/mV2) 
  d2l = -3.017916391902783e-08     (/mV3) 

  an = 0.05795847284739256     (/mV) 
  bn = 0.6520582761660094     (1) 
  vhn = -42.301211060330765     (mV) 
  An = 0.7513078725176383     (/ms) 
  b1n = 0.09808502009804915     (/mV) 
  c1n = 0.0011532937076404724     (/mV2) 
  d1n = -1.1952436753184095e-05     (/mV3) 
  b2n = -0.022261546191178004     (/mV) 
  c2n = 0.00025257848072110766     (/mV2) 
  d2n = -3.968523792694852e-08     (/mV3) 
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