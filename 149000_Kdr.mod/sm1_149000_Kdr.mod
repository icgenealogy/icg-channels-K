NEURON
{
  SUFFIX Kdr 
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

  ah = -0.24919138133338523     (/mV) 
  bh = 6.2326247047306795     (1) 
  vhh = -21.153146251453286     (mV) 
  Ah = 1200.2216140270357     (/ms) 
  b1h = 1.4299502913738298e-05     (/mV) 
  c1h = 3.25600817093369e-07     (/mV2) 
  d1h = 2.2201798242299073e-09     (/mV3) 
  b2h = 4.159193913284825     (/mV) 
  c2h = -0.013451191287801716     (/mV2) 
  d2h = 1.455244390581005e-05     (/mV3) 

  am = 0.0848676286449557     (/mV) 
  bm = -0.9853973876876436     (1) 
  vhm = -31.072642200541548     (mV) 
  Am = 45.49513071744733     (/ms) 
  b1m = -0.05558187056303903     (/mV) 
  c1m = 0.0003775715184894439     (/mV2) 
  d1m = -1.0621446380659863e-06     (/mV3) 
  b2m = -0.03259266711235619     (/mV) 
  c2m = -1.238990599988242e-05     (/mV2) 
  d2m = -2.0013806426943047e-07     (/mV3) 
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
  g = gbar*h*m*m
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