NEURON
{
  SUFFIX kf 
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

  an = 0.03448271053806348     (/mV) 
  bn = -1.6206865226498504     (1) 
  vhn = 20.869903730627374     (mV) 
  An = 0.2523617371731788     (/ms) 
  b1n = -0.010169300921854568     (/mV) 
  c1n = -2.8859212209199275e-05     (/mV2) 
  d1n = 3.342092600042666e-07     (/mV3) 
  b2n = 0.010169368627728667     (/mV) 
  c2n = -0.0003002199836218994     (/mV2) 
  d2n = -2.665789528809524e-07     (/mV3) 

  al = -0.09999882249184666     (/mV) 
  bl = 6.599940463663558     (1) 
  vhl = -70.74636039448148     (mV) 
  Al = 32.79905204106816     (/ms) 
  b1l = 0.07751033853445996     (/mV) 
  c1l = 0.0020211382554365647     (/mV2) 
  d1l = 4.849468601595324e-05     (/mV3) 
  b2l = 0.08867326261821445     (/mV) 
  c2l = -0.0008809113136608224     (/mV2) 
  d2l = 2.6961653262206927e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
  lInf 
  lTau 
}

STATE
{
  n
  l
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*l
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
  l' = (lInf - l) / lTau 
}

INITIAL
{
  rates(v)
  n = nInf 
  l = lInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 


  UNITSON
}