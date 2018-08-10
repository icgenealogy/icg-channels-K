NEURON
{
  SUFFIX iA 
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

  am1 = 0.11764693232294225     (/mV) 
  bm1 = -7.0588156643147135     (1) 
  vhm1 = -94.19158854450902     (mV) 
  Am1 = 0.14347999015693644     (/ms) 
  b1m1 = 0.05945506780682384     (/mV) 
  c1m1 = 0.0009153997631916628     (/mV2) 
  d1m1 = -1.2247158930157668e-05     (/mV3) 
  b2m1 = -0.059452950868179984     (/mV) 
  c2m1 = 0.0012868251283380966     (/mV2) 
  d2m1 = -5.434391457018805e-06     (/mV3) 

  ah1 = -0.1666645706056242     (/mV) 
  bh1 = 12.99985888133849     (1) 
  vhh1 = -69.99641598269255     (mV) 
  Ah1 = 21.16283989972838     (/ms) 
  b1h1 = -0.09750809763864891     (/mV) 
  c1h1 = 0.0016976836957911258     (/mV2) 
  d1h1 = -6.883440859195506e-06     (/mV3) 
  b2h1 = -0.012659697260573746     (/mV) 
  c2h1 = 0.0005070296515200169     (/mV2) 
  d2h1 = -2.409052794096825e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  m1Inf 
  m1Tau 
  h1Inf 
  h1Tau 
}

STATE
{
  m1
  h1
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m1*h1
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  m1' = (m1Inf - m1) / m1Tau 
  h1' = (h1Inf - h1) / h1Tau 
}

INITIAL
{
  rates(v)
  m1 = m1Inf 
  h1 = h1Inf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    m1Inf = 1/(1 + exp(-am1*v + bm1)) 
    m1Tau = Am1 / ( exp(-(b1m1*(v-vhm1) + c1m1*(v-vhm1)^2 + d1m1*(v-vhm1)^3)) + exp((b2m1*(v-vhm1) + c2m1*(v-vhm1)^2 + d2m1*(v-vhm1)^3)) ) 

    h1Inf = 1/(1 + exp(-ah1*v + bh1)) 
    h1Tau = Ah1 / ( exp(-(b1h1*(v-vhh1) + c1h1*(v-vhh1)^2 + d1h1*(v-vhh1)^3)) + exp((b2h1*(v-vhh1) + c2h1*(v-vhh1)^2 + d2h1*(v-vhh1)^3)) ) 


  UNITSON
}