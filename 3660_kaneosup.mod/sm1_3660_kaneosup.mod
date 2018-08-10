NEURON
{
  SUFFIX kasup 
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

  ah = -0.1265818829226449     (/mV) 
  bh = 6.493655648132147     (1) 
  vhh = -68.91482285500334     (mV) 
  Ah = 11.948137673421531     (/ms) 
  b1h = 0.17141566000302624     (/mV) 
  c1h = -0.002045952408753418     (/mV2) 
  d1h = 6.919637480965495e-06     (/mV3) 
  b2h = 0.006069524957867182     (/mV) 
  c2h = -4.79628365188711e-05     (/mV2) 
  d2h = 1.6245472875509395e-07     (/mV3) 

  am = 0.0862067680268906     (/mV) 
  bm = -0.0689563775546829     (1) 
  vhm = -44.11799383947714     (mV) 
  Am = 4.926333066177332     (/ms) 
  b1m = 0.07148962602659088     (/mV) 
  c1m = -0.0001363970315708123     (/mV2) 
  d1m = -2.4819839254868154e-06     (/mV3) 
  b2m = 0.02608556373886577     (/mV) 
  c2m = -8.840413929814383e-06     (/mV2) 
  d2m = 1.0064839378278515e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}