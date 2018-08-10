NEURON
{
  SUFFIX IA 
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

  aa = 0.06024084499003787     (/mV) 
  ba = -0.8433673713459108     (1) 
  vha = -73.56183466673515     (mV) 
  Aa = 10.020000007078417     (/ms) 
  b1a = -3.322440488619548e-09     (/mV) 
  c1a = 1.803021258416065e-09     (/mV2) 
  d1a = -5.310905572541754e-11     (/mV3) 
  b2a = -3.3883273721517777e-09     (/mV) 
  c2a = 1.8063835814949279e-09     (/mV2) 
  d2a = -5.3134583518005927e-11     (/mV3) 

  ab = -0.13694745381539478     (/mV) 
  bb = 9.723552197742857     (1) 
  vhb = -88.450423433765     (mV) 
  Ab = 288.8013970042707     (/ms) 
  b1b = -0.08708826344573797     (/mV) 
  c1b = 0.0007668362511210069     (/mV2) 
  d1b = -2.100146748661933e-06     (/mV3) 
  b2b = -0.09413528756393619     (/mV) 
  c2b = -0.003615377538042881     (/mV2) 
  d2b = -0.00014631928810300026     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  aInf 
  aTau 
  bInf 
  bTau 
}

STATE
{
  a
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*a*b
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  a' = (aInf - a) / aTau 
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  a = aInf 
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    aInf = 1/(1 + exp(-aa*v + ba)) 
    aTau = Aa / ( exp(-(b1a*(v-vha) + c1a*(v-vha)^2 + d1a*(v-vha)^3)) + exp((b2a*(v-vha) + c2a*(v-vha)^2 + d2a*(v-vha)^3)) ) 

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}