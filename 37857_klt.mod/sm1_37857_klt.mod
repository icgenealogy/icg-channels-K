NEURON
{
  SUFFIX klt 
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

  aw = 0.08046235370827627     (/mV) 
  bw = -5.35666536606703     (1) 
  vhw = -70.10209966783464     (mV) 
  Aw = 2.4680627930291132     (/ms) 
  b1w = -0.07485762660815536     (/mV) 
  c1w = 0.0007668683133353895     (/mV2) 
  d1w = -2.4000171198097527e-06     (/mV3) 
  b2w = -0.10037814143155624     (/mV) 
  c2w = -0.004009023492916943     (/mV2) 
  d2w = -6.455584995212921e-05     (/mV3) 

  az = -0.00719279408515454     (/mV) 
  bz = -0.31284387836803373     (1) 
  vhz = -59.3058377503388     (mV) 
  Az = 216.14697164901565     (/ms) 
  b1z = -0.05763659677637189     (/mV) 
  c1z = 0.00032309529165728835     (/mV2) 
  d1z = -4.939896935705804e-07     (/mV3) 
  b2z = -0.12100415663917173     (/mV) 
  c2z = -0.0007500508267301478     (/mV2) 
  d2z = 9.749623431472195e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  wInf 
  wTau 
  zInf 
  zTau 
}

STATE
{
  w
  z
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*w*w*w*w*z
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  w' = (wInf - w) / wTau 
  z' = (zInf - z) / zTau 
}

INITIAL
{
  rates(v)
  w = wInf 
  z = zInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    wInf = 1/(1 + exp(-aw*v + bw)) 
    wTau = Aw / ( exp(-(b1w*(v-vhw) + c1w*(v-vhw)^2 + d1w*(v-vhw)^3)) + exp((b2w*(v-vhw) + c2w*(v-vhw)^2 + d2w*(v-vhw)^3)) ) 

    zInf = 1/(1 + exp(-az*v + bz)) 
    zTau = Az / ( exp(-(b1z*(v-vhz) + c1z*(v-vhz)^2 + d1z*(v-vhz)^3)) + exp((b2z*(v-vhz) + c2z*(v-vhz)^2 + d2z*(v-vhz)^3)) ) 


  UNITSON
}