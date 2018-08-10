NEURON
{
  SUFFIX kfast 
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

  aa = 0.03845920112430915     (/mV) 
  ba = -2.15371871643928     (1) 
  vha = -41.1001218021942     (mV) 
  Aa = 0.6011581360607032     (/ms) 
  b1a = 0.17146759662049305     (/mV) 
  c1a = 0.0052839644515699025     (/mV2) 
  d1a = 4.94791089952615e-05     (/mV3) 
  b2a = 0.025020095793627354     (/mV) 
  c2a = -0.00016612893552304253     (/mV2) 
  d2a = 4.859154722860598e-07     (/mV3) 

  ab = -0.09086578982450015     (/mV) 
  bb = 6.815082117471579     (1) 
  vhb = -83.98211149791268     (mV) 
  Ab = 210.34802542813566     (/ms) 
  b1b = -0.1515422303368159     (/mV) 
  c1b = 0.0016604169131034208     (/mV2) 
  d1b = -5.291776055862826e-06     (/mV3) 
  b2b = -0.2779626336120973     (/mV) 
  c2b = -0.025634031963282067     (/mV2) 
  d2b = -0.0008397465340210808     (/mV3) 
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