NEURON
{
  SUFFIX GrC_KA 
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

  aa = 0.05050496183659363     (/mV) 
  ba = -2.358579655594201     (1) 
  vha = -89.99990836734331     (mV) 
  Aa = 7.247177591160779     (/ms) 
  b1a = -70.57672665551853     (/mV) 
  c1a = 0.005404858571671471     (/mV2) 
  d1a = -1.7097226776267418e-05     (/mV3) 
  b2a = -69.50969433604297     (/mV) 
  c2a = -0.0057495519864227155     (/mV2) 
  d2a = 1.68806865809225e-05     (/mV3) 

  ab = -0.11904711458907849     (/mV) 
  bb = 9.380924925052406     (1) 
  vhb = -79.36854555455744     (mV) 
  Ab = 25.68719854142279     (/ms) 
  b1b = 0.13580722485449787     (/mV) 
  c1b = 0.004864496338589226     (/mV2) 
  d1b = 0.00010890758463996085     (/mV3) 
  b2b = 0.09380761800292697     (/mV) 
  c2b = -0.0009111443946453374     (/mV2) 
  d2b = 2.708777900264533e-06     (/mV3) 
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