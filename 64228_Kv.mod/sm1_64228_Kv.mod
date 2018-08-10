NEURON
{
  SUFFIX Kv 
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

  amKv = 0.038099520091651566     (/mV) 
  bmKv = -0.00789305137788822     (1) 
  vhmKv = -41.42251486593902     (mV) 
  AmKv = 16.031414328230017     (/ms) 
  b1mKv = 0.015416775701911068     (/mV) 
  c1mKv = 8.253803776160178e-06     (/mV2) 
  d1mKv = -1.021054315440223e-06     (/mV3) 
  b2mKv = 0.004711122443185555     (/mV) 
  c2mKv = 0.00016597862449561045     (/mV2) 
  d2mKv = -8.283014900758373e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mKvInf 
  mKvTau 
}

STATE
{
  mKv
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*mKv*mKv*mKv*mKv
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  mKv' = (mKvInf - mKv) / mKvTau 
}

INITIAL
{
  rates(v)
  mKv = mKvInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mKvInf = 1/(1 + exp(-amKv*v + bmKv)) 
    mKvTau = AmKv / ( exp(-(b1mKv*(v-vhmKv) + c1mKv*(v-vhmKv)^2 + d1mKv*(v-vhmKv)^3)) + exp((b2mKv*(v-vhmKv) + c2mKv*(v-vhmKv)^2 + d2mKv*(v-vhmKv)^3)) ) 


  UNITSON
}