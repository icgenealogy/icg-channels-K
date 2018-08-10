NEURON
{
  SUFFIX IKaSM 
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

  ah = -0.09999877589344243     (/mV) 
  bh = 8.599922261162419     (1) 
  vhh = -251.03320368888714     (mV) 
  Ah = 44.99123417841185     (/ms) 
  b1h = -0.005298804794553143     (/mV) 
  c1h = 1.2349928688627149e-05     (/mV2) 
  d1h = -1.0097219161084487e-08     (/mV3) 
  b2h = -0.0006243446188990025     (/mV) 
  c2h = -3.0393455691266166e-05     (/mV2) 
  d2h = 3.2638868581027335e-08     (/mV3) 

  am = 0.10223436702266264     (/mV) 
  bm = -5.4869731812862605     (1) 
  vhm = -30.27971714662122     (mV) 
  Am = 1.6512588961078039     (/ms) 
  b1m = -0.00792272132886791     (/mV) 
  c1m = -0.000286125371857165     (/mV2) 
  d1m = 2.013648817209463e-06     (/mV3) 
  b2m = 0.039086553691600544     (/mV) 
  c2m = -0.0007042379474484638     (/mV2) 
  d2m = 3.768387457582591e-06     (/mV3) 
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
  g = gbar*h*m*m*m*m
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