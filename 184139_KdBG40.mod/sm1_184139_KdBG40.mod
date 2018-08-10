NEURON
{
  SUFFIX KdBG 
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

  axs = 0.19499981176164427     (/mV) 
  bxs = -9.359985615628519     (1) 
  vhxs = -5.621991789615483     (mV) 
  Axs = 2.1259375851416253     (/ms) 
  b1xs = -0.005161522290973455     (/mV) 
  c1xs = -7.383638727223121e-05     (/mV2) 
  d1xs = 5.162721739759551e-07     (/mV3) 
  b2xs = 0.005233143924347421     (/mV) 
  c2xs = -0.00017087543393690783     (/mV2) 
  d2xs = -1.8754537646513987e-06     (/mV3) 

  ays = -0.1169406765927602     (/mV) 
  bys = 10.525959592604746     (1) 
  vhys = -78.92016519297218     (mV) 
  Ays = 1018.615417393228     (/ms) 
  b1ys = 0.016420553530714505     (/mV) 
  c1ys = -0.0004807912237656121     (/mV2) 
  d1ys = 1.981205709126352e-06     (/mV3) 
  b2ys = 0.06989524732802463     (/mV) 
  c2ys = -0.0007409680027823032     (/mV2) 
  d2ys = 2.30747186941109e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  xsInf 
  xsTau 
  ysInf 
  ysTau 
}

STATE
{
  xs
  ys
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*xs*ys
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  xs' = (xsInf - xs) / xsTau 
  ys' = (ysInf - ys) / ysTau 
}

INITIAL
{
  rates(v)
  xs = xsInf 
  ys = ysInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    xsInf = 1/(1 + exp(-axs*v + bxs)) 
    xsTau = Axs / ( exp(-(b1xs*(v-vhxs) + c1xs*(v-vhxs)^2 + d1xs*(v-vhxs)^3)) + exp((b2xs*(v-vhxs) + c2xs*(v-vhxs)^2 + d2xs*(v-vhxs)^3)) ) 

    ysInf = 1/(1 + exp(-ays*v + bys)) 
    ysTau = Ays / ( exp(-(b1ys*(v-vhys) + c1ys*(v-vhys)^2 + d1ys*(v-vhys)^3)) + exp((b2ys*(v-vhys) + c2ys*(v-vhys)^2 + d2ys*(v-vhys)^3)) ) 


  UNITSON
}