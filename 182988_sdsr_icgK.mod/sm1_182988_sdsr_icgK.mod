NEURON
{
  SUFFIX sdsr 
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

  aasd = 0.09999973126596025     (/mV) 
  basd = -3.9999820583267276     (1) 
  vhasd = -46.69413497772626     (mV) 
  Aasd = 5.360000002236288     (/ms) 
  b1asd = 1.7694449145534415e-08     (/mV) 
  c1asd = 4.070021865019931e-09     (/mV2) 
  d1asd = 3.239967056968791e-11     (/mV3) 
  b2asd = 1.7734839958276396e-08     (/mV) 
  c2asd = 4.070860013811846e-09     (/mV2) 
  d2asd = 3.2381573228246595e-11     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  asdInf 
  asdTau 
}

STATE
{
  asd
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*asd
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  asd' = (asdInf - asd) / asdTau 
}

INITIAL
{
  rates(v)
  asd = asdInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    asdInf = 1/(1 + exp(-aasd*v + basd)) 
    asdTau = Aasd / ( exp(-(b1asd*(v-vhasd) + c1asd*(v-vhasd)^2 + d1asd*(v-vhasd)^3)) + exp((b2asd*(v-vhasd) + c2asd*(v-vhasd)^2 + d2asd*(v-vhasd)^3)) ) 


  UNITSON
}