NEURON
{
  SUFFIX Ikdrs 
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

  ah = -0.053409673038965544     (/mV) 
  bh = 2.588708673463563     (1) 
  vhh = -80.62936232610028     (mV) 
  Ah = 2000.0120179992857     (/ms) 
  b1h = 0.0028452666291124825     (/mV) 
  c1h = 3.890905846241037e-06     (/mV2) 
  d1h = 1.2777582948202538e-09     (/mV3) 
  b2h = 0.0028457705410125404     (/mV) 
  c2h = -4.215000460935736e-06     (/mV2) 
  d2h = 2.304946847040307e-09     (/mV3) 

  am = 0.056804919856396144     (/mV) 
  bm = -0.0966719825882859     (1) 
  vhm = 8.304461967603432     (mV) 
  Am = 14.010018528698991     (/ms) 
  b1m = 0.018546279183489264     (/mV) 
  c1m = 0.0001705651512694397     (/mV2) 
  d1m = 5.360025137301647e-07     (/mV3) 
  b2m = 0.01872202437629011     (/mV) 
  c2m = -0.00017824422015373935     (/mV2) 
  d2m = 5.996453123513909e-07     (/mV3) 
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