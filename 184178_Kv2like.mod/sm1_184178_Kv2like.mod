NEURON
{
  SUFFIX Kv2like 
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

  ah1 = -0.0909089494286563     (/mV) 
  bh1 = 5.2727178280249705     (1) 
  vhh1 = -17.387434397008143     (mV) 
  Ah1 = 423.57505168823934     (/ms) 
  b1h1 = -0.008701232477046154     (/mV) 
  c1h1 = -0.00021728522520302686     (/mV2) 
  d1h1 = 1.572724035856965e-06     (/mV3) 
  b2h1 = 0.03161509392018863     (/mV) 
  c2h1 = -0.00017740345348861024     (/mV2) 
  d2h1 = -3.1861512753814757e-06     (/mV3) 

  am = 0.08299548012487284     (/mV) 
  bm = -1.69260260840302     (1) 
  vhm = -18.243811287866095     (mV) 
  Am = 25.75616601252719     (/ms) 
  b1m = 0.004126740997746085     (/mV) 
  c1m = -5.777704624096709e-05     (/mV2) 
  d1m = -2.702254136301662e-07     (/mV3) 
  b2m = 0.079238434408131     (/mV) 
  c2m = -0.00020503427096014572     (/mV2) 
  d2m = -6.669138816171254e-07     (/mV3) 

  ah2 = -0.09076493561251471     (/mV) 
  bh2 = 5.267046434379604     (1) 
  vhh2 = -102.11770477887438     (mV) 
  Ah2 = 1115.0141466825257     (/ms) 
  b1h2 = -0.0005380885086383229     (/mV) 
  c1h2 = -1.7326282147064914e-05     (/mV2) 
  d1h2 = 4.90747603119953e-08     (/mV3) 
  b2h2 = -0.07697526549618622     (/mV) 
  c2h2 = 0.0006928752029651623     (/mV2) 
  d2h2 = -7.374469824919328e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  h1Inf 
  h1Tau 
  mInf 
  mTau 
  h2Inf 
  h2Tau 
}

STATE
{
  h1
  m
  h2
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h1*m*m*h2
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h1' = (h1Inf - h1) / h1Tau 
  m' = (mInf - m) / mTau 
  h2' = (h2Inf - h2) / h2Tau 
}

INITIAL
{
  rates(v)
  h1 = h1Inf 
  m = mInf 
  h2 = h2Inf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    h1Inf = 1/(1 + exp(-ah1*v + bh1)) 
    h1Tau = Ah1 / ( exp(-(b1h1*(v-vhh1) + c1h1*(v-vhh1)^2 + d1h1*(v-vhh1)^3)) + exp((b2h1*(v-vhh1) + c2h1*(v-vhh1)^2 + d2h1*(v-vhh1)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    h2Inf = 1/(1 + exp(-ah2*v + bh2)) 
    h2Tau = Ah2 / ( exp(-(b1h2*(v-vhh2) + c1h2*(v-vhh2)^2 + d1h2*(v-vhh2)^3)) + exp((b2h2*(v-vhh2) + c2h2*(v-vhh2)^2 + d2h2*(v-vhh2)^3)) ) 


  UNITSON
}