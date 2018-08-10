NEURON
{
  SUFFIX GRC_KA 
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

  aa = 0.05858555593985653     (/mV) 
  ba = -2.2253278557895153     (1) 
  vha = -23.337190963029002     (mV) 
  Aa = 1.14977024572074     (/ms) 
  b1a = 0.0468768905083587     (/mV) 
  c1a = -2.443035689722485e-05     (/mV2) 
  d1a = -5.496344306714484e-07     (/mV3) 
  b2a = 0.03475545432774884     (/mV) 
  c2a = -0.0003846241063304799     (/mV2) 
  d2a = 1.3876136422559905e-06     (/mV3) 

  ab = -0.11904233611883472     (/mV) 
  bb = 9.380661094424406     (1) 
  vhb = -79.33111316816574     (mV) 
  Ab = 141.02777789783036     (/ms) 
  b1b = 0.13543268762832722     (/mV) 
  c1b = 0.004833285334705666     (/mV2) 
  d1b = 0.00010826006803148711     (/mV3) 
  b2b = 0.0940390754962247     (/mV) 
  c2b = -0.0009139223744927708     (/mV2) 
  d2b = 2.7182016339962352e-06     (/mV3) 
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
  g = gbar*a*a*a*b
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