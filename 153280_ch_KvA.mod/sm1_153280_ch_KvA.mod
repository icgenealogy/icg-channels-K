NEURON
{
  SUFFIX ch_KvA 
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

  an = 0.11223026645866716     (/mV) 
  bn = -3.7709400395677646     (1) 
  vhn = -33.734694997607455     (mV) 
  An = 23.21567306248117     (/ms) 
  b1n = -0.06671708866425222     (/mV) 
  c1n = -1.8797195281082652e-05     (/mV2) 
  d1n = 2.3822693226601747e-07     (/mV3) 
  b2n = -0.04514696562842057     (/mV) 
  c2n = -1.808223729488769e-06     (/mV2) 
  d2n = 1.987126227915805e-08     (/mV3) 

  al = -0.14964023030337495     (/mV) 
  bl = 12.420143539444984     (1) 
  vhl = -77.35463531193386     (mV) 
  Al = 8.069252434639186     (/ms) 
  b1l = 0.09247447159362802     (/mV) 
  c1l = -0.0016681773302173423     (/mV2) 
  d1l = 1.4852669527109243e-05     (/mV3) 
  b2l = 0.007847077690378294     (/mV) 
  c2l = -6.0532062138894286e-05     (/mV2) 
  d2l = 1.518256902902515e-07     (/mV3) 
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