NEURON
{
  SUFFIX ka 
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

  aa = 0.07806735972650977     (/mV) 
  ba = -3.8855187579758197     (1) 
  vha = -41.93887458777563     (mV) 
  Aa = 1.0327381416943526     (/ms) 
  b1a = -0.08921385645477775     (/mV) 
  c1a = 0.0006326336631683502     (/mV2) 
  d1a = -1.260289705714064e-06     (/mV3) 
  b2a = -0.025958136706417625     (/mV) 
  c2a = 0.00019018371560636918     (/mV2) 
  d2a = 1.0511889612799972e-06     (/mV3) 

  ac = -0.10517586334518683     (/mV) 
  bc = 6.039175967479708     (1) 
  vhc = 48.33735186771688     (mV) 
  Ac = 38.422322541683215     (/ms) 
  b1c = 0.0002494434710755949     (/mV) 
  c1c = -4.388994034797139e-05     (/mV2) 
  d1c = 4.104706067973172e-07     (/mV3) 
  b2c = -0.0002641976029176169     (/mV) 
  c2c = -4.355590280643075e-05     (/mV2) 
  d2c = 6.096308532455468e-07     (/mV3) 

  ab = -0.10517603033003764     (/mV) 
  bb = 6.039173943753242     (1) 
  vhb = -48.678206865951815     (mV) 
  Ab = 10.162172017406176     (/ms) 
  b1b = 0.03622186390207691     (/mV) 
  c1b = -0.00010507670697075406     (/mV2) 
  d1b = -1.4651284525811542e-06     (/mV3) 
  b2b = 0.03786757825773754     (/mV) 
  c2b = -3.3960402826648556e-06     (/mV2) 
  d2b = -5.792524940743331e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  aInf 
  aTau 
  cInf 
  cTau 
  bInf 
  bTau 
}

STATE
{
  a
  c
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*a*a*a*a*c*b
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  a' = (aInf - a) / aTau 
  c' = (cInf - c) / cTau 
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  a = aInf 
  c = cInf 
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    aInf = 1/(1 + exp(-aa*v + ba)) 
    aTau = Aa / ( exp(-(b1a*(v-vha) + c1a*(v-vha)^2 + d1a*(v-vha)^3)) + exp((b2a*(v-vha) + c2a*(v-vha)^2 + d2a*(v-vha)^3)) ) 

    cInf = 1/(1 + exp(-ac*v + bc)) 
    cTau = Ac / ( exp(-(b1c*(v-vhc) + c1c*(v-vhc)^2 + d1c*(v-vhc)^3)) + exp((b2c*(v-vhc) + c2c*(v-vhc)^2 + d2c*(v-vhc)^3)) ) 

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}