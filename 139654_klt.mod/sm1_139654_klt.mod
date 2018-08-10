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

  aw = 0.08041137549842872     (/mV) 
  bw = -5.352810712530741     (1) 
  vhw = -70.1026465364745     (mV) 
  Aw = 6.3822546657067605     (/ms) 
  b1w = -0.0749634751038156     (/mV) 
  c1w = 0.0007708709501078465     (/mV2) 
  d1w = -2.417836866977962e-06     (/mV3) 
  b2w = -0.10073346737336143     (/mV) 
  c2w = -0.0040557631995864285     (/mV2) 
  d2w = -6.585733877442789e-05     (/mV3) 
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
}

STATE
{
  w
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*w*w*w*w
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  w' = (wInf - w) / wTau 
}

INITIAL
{
  rates(v)
  w = wInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    wInf = 1/(1 + exp(-aw*v + bw)) 
    wTau = Aw / ( exp(-(b1w*(v-vhw) + c1w*(v-vhw)^2 + d1w*(v-vhw)^3)) + exp((b2w*(v-vhw) + c2w*(v-vhw)^2 + d2w*(v-vhw)^3)) ) 


  UNITSON
}