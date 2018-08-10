NEURON
{
  SUFFIX kad 
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

  an = 0.07019377066189549     (/mV) 
  bn = -0.0495089366400468     (1) 
  vhn = -23.303552532833553     (mV) 
  An = 0.7659245108493328     (/ms) 
  b1n = 0.0028188187070493204     (/mV) 
  c1n = -0.00024317938788117916     (/mV2) 
  d1n = 6.292860337563349e-07     (/mV3) 
  b2n = -0.07851481636260409     (/mV) 
  c2n = -0.00046689393897065266     (/mV2) 
  d2n = 7.961544273526441e-06     (/mV3) 

  al = -0.11223025470131971     (/mV) 
  bl = 6.2848930203926265     (1) 
  vhl = -35.88305559744569     (mV) 
  Al = 7.427560117926344     (/ms) 
  b1l = 0.07567915523058798     (/mV) 
  c1l = 0.0012991466120830881     (/mV2) 
  d1l = -1.0210684200064545e-05     (/mV3) 
  b2l = -0.01898804379987625     (/mV) 
  c2l = 5.422340271537557e-05     (/mV2) 
  d2l = -3.014454781007447e-08     (/mV3) 
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