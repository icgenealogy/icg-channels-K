NEURON
{
  SUFFIX iA 
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

  an = 0.03448268184397103     (/mV) 
  bn = -1.620677338177199     (1) 
  vhn = 21.779496247144724     (mV) 
  An = 0.8433072855419395     (/ms) 
  b1n = -0.010306632913526736     (/mV) 
  c1n = 0.0003186175649428213     (/mV2) 
  d1n = 5.028550579954935e-07     (/mV3) 
  b2n = 0.010306242188336796     (/mV) 
  c2n = 2.3809494751467786e-05     (/mV2) 
  d2n = -3.6462822360262446e-07     (/mV3) 

  al = -0.09999199695196338     (/mV) 
  bl = 6.599563388614644     (1) 
  vhl = -70.7298181510103     (mV) 
  Al = 114.36884118293675     (/ms) 
  b1l = -0.08872347024976067     (/mV) 
  c1l = 0.0008814857615237269     (/mV2) 
  d1l = -2.698151679597354e-06     (/mV3) 
  b2l = -0.07739857608842408     (/mV) 
  c2l = -0.002015673552090948     (/mV2) 
  d2l = -4.8407633775241495e-05     (/mV3) 
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