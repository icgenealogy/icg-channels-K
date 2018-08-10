NEURON
{
  SUFFIX a 
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

  aa = 0.17499967391667295     (/mV) 
  ba = -11.374978805205194     (1) 
  vha = -63.33130924636133     (mV) 
  Aa = 4.337018515493505     (/ms) 
  b1a = 0.011015889385961356     (/mV) 
  c1a = -0.0004833955800735504     (/mV2) 
  d1a = 2.4070853878675598e-06     (/mV3) 
  b2a = 0.09015400988586549     (/mV) 
  c2a = -0.0016658885426057244     (/mV2) 
  d2a = 7.077013485304588e-06     (/mV3) 

  ab = -0.27398236059014014     (/mV) 
  bb = 19.452852180136517     (1) 
  vhb = -58.34803985843144     (mV) 
  Ab = 61.75129838427687     (/ms) 
  b1b = -0.049941717582043535     (/mV) 
  c1b = 0.0011144168712679424     (/mV2) 
  d1b = -5.148099347577716e-06     (/mV3) 
  b2b = -0.010098509840965546     (/mV) 
  c2b = 0.00037033475946549017     (/mV2) 
  d2b = -1.8770047260062435e-06     (/mV3) 
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
  bInf 
  bTau 
}

STATE
{
  a
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*a*b
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  a' = (aInf - a) / aTau 
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  a = aInf 
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    aInf = 1/(1 + exp(-aa*v + ba)) 
    aTau = Aa / ( exp(-(b1a*(v-vha) + c1a*(v-vha)^2 + d1a*(v-vha)^3)) + exp((b2a*(v-vha) + c2a*(v-vha)^2 + d2a*(v-vha)^3)) ) 

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}