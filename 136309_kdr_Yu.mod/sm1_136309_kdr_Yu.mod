NEURON
{
  SUFFIX kdrYu 
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

  am = 0.11024079268095185     (/mV) 
  bm = 0.47989224513772605     (1) 
  vhm = 11.257972663612954     (mV) 
  Am = 17.888013465500922     (/ms) 
  b1m = 0.009631821823572     (/mV) 
  c1m = -0.00021783886865238023     (/mV2) 
  d1m = -1.8262134843722947e-06     (/mV3) 
  b2m = 0.06899265485431579     (/mV) 
  c2m = -0.0002208298427994608     (/mV2) 
  d2m = -7.576033725241335e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}