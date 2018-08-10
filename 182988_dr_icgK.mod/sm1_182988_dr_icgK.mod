NEURON
{
  SUFFIX dr 
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

  aad = 0.24999489569329322     (/mV) 
  bad = -6.249814903408923     (1) 
  vhad = -158.39321243553206     (mV) 
  Aad = 0.4436009711336736     (/ms) 
  b1ad = -0.002377693815465407     (/mV) 
  c1ad = 2.245109241377215e-05     (/mV2) 
  d1ad = -6.675379936157486e-08     (/mV3) 
  b2ad = -1.220664345863341     (/mV) 
  c2ad = -0.0071187801622626485     (/mV2) 
  d2ad = 6.486912676576076e-05     (/mV3) 

  aar = 0.24999489569328584     (/mV) 
  bar = -6.249814903408743     (1) 
  vhar = -64.50921006453557     (mV) 
  Aar = 0.82000000031825     (/ms) 
  b1ar = -2.6824613215787995e-09     (/mV) 
  c1ar = -7.479794048618068e-10     (/mV2) 
  d1ar = 3.6325597393886215e-11     (/mV3) 
  b2ar = -2.6882558735958286e-09     (/mV) 
  c2ar = -7.468506516935005e-10     (/mV2) 
  d2ar = 3.631513343684531e-11     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  adInf 
  adTau 
  arInf 
  arTau 
}

STATE
{
  ad
  ar
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*ar
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  ad' = (adInf - ad) / adTau 
  ar' = (arInf - ar) / arTau 
}

INITIAL
{
  rates(v)
  ad = adInf 
  ar = arInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    adInf = 1/(1 + exp(-aad*v + bad)) 
    adTau = Aad / ( exp(-(b1ad*(v-vhad) + c1ad*(v-vhad)^2 + d1ad*(v-vhad)^3)) + exp((b2ad*(v-vhad) + c2ad*(v-vhad)^2 + d2ad*(v-vhad)^3)) ) 

    arInf = 1/(1 + exp(-aar*v + bar)) 
    arTau = Aar / ( exp(-(b1ar*(v-vhar) + c1ar*(v-vhar)^2 + d1ar*(v-vhar)^3)) + exp((b2ar*(v-vhar) + c2ar*(v-vhar)^2 + d2ar*(v-vhar)^3)) ) 


  UNITSON
}