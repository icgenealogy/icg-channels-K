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

  al = -0.11223019851510646     (/mV) 
  bl = 5.162585389585405     (1) 
  vhl = 143.76510376452111     (mV) 
  Al = 42.69645920838157     (/ms) 
  b1l = -0.3150353488488659     (/mV) 
  c1l = -0.003455062532851364     (/mV2) 
  d1l = -9.061825204826743e-06     (/mV3) 
  b2l = -0.0011904365784350546     (/mV) 
  c2l = 6.028165606564007e-05     (/mV2) 
  d2l = 4.378404860562233e-08     (/mV3) 

  an = 0.057964773992890364     (/mV) 
  bn = 1.2317494646315688     (1) 
  vhn = -30.875531024463736     (mV) 
  An = 0.8243818413222268     (/ms) 
  b1n = 0.018656380544677723     (/mV) 
  c1n = -0.00020509713452935237     (/mV2) 
  d1n = -2.1467675860825926e-07     (/mV3) 
  b2n = -0.09716251070732923     (/mV) 
  c2n = -0.000910011055752918     (/mV2) 
  d2n = 1.1050514606090258e-05     (/mV3) 
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