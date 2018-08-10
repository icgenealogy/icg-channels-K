NEURON
{
  SUFFIX kvR213W 
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

  am = 0.06565804662216644     (/mV) 
  bm = -1.454309817395866     (1) 
  vhm = -10.638760657498844     (mV) 
  Am = 23.08364343080579     (/ms) 
  b1m = 0.0005484805806574222     (/mV) 
  c1m = -0.0001758612952624846     (/mV2) 
  d1m = 1.3717319827110862e-06     (/mV3) 
  b2m = -0.08436160695139358     (/mV) 
  c2m = -0.0011877817193029893     (/mV2) 
  d2m = 1.622710467489937e-05     (/mV3) 
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