NEURON
{
  SUFFIX GRC_KIR 
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

  ad = -0.06899986349282751     (/mV) 
  bd = 6.037774761399099     (1) 
  vhd = -22.920025562114432     (mV) 
  Ad = 0.3439898463049937     (/ms) 
  b1d = 0.01052873925216924     (/mV) 
  c1d = 0.0007089500384939453     (/mV2) 
  d1d = 6.9588736897411345e-06     (/mV3) 
  b2d = 0.06295039262973778     (/mV) 
  c2d = -0.0006525213671411371     (/mV2) 
  d2d = 2.1723041817384198e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  dInf 
  dTau 
}

STATE
{
  d
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*d
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  d' = (dInf - d) / dTau 
}

INITIAL
{
  rates(v)
  d = dInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    dInf = 1/(1 + exp(-ad*v + bd)) 
    dTau = Ad / ( exp(-(b1d*(v-vhd) + c1d*(v-vhd)^2 + d1d*(v-vhd)^3)) + exp((b2d*(v-vhd) + c2d*(v-vhd)^2 + d2d*(v-vhd)^3)) ) 


  UNITSON
}