NEURON
{
  SUFFIX spike 
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

  an = 0.050937544942293606     (/mV) 
  bn = -1.3822129653088513     (1) 
  vhn = -62.21513700641488     (mV) 
  An = 3.8534045899263836     (/ms) 
  b1n = -0.029067099260923154     (/mV) 
  c1n = 0.00012934272405133883     (/mV2) 
  d1n = -2.5312683267826506e-07     (/mV3) 
  b2n = -0.037665835111076434     (/mV) 
  c2n = -0.0004572783414851643     (/mV2) 
  d2n = -3.466498634724396e-06     (/mV3) 

  ah = -0.14173676245080263     (/mV) 
  bh = 6.732956795375899     (1) 
  vhh = -50.090019440001896     (mV) 
  Ah = 2.9599163006516713     (/ms) 
  b1h = 0.07536146717436855     (/mV) 
  c1h = 0.0010494398553990507     (/mV2) 
  d1h = 1.2552233057621863e-05     (/mV3) 
  b2h = 0.09956799624219982     (/mV) 
  c2h = -0.0010698658570481665     (/mV2) 
  d2h = 3.6034284914899084e-06     (/mV3) 

  ap = 0.1199919006249093     (/mV) 
  bp = -4.890074438654973     (1) 
  vhp = -21.608292284396022     (mV) 
  Ap = 4.43366996595335     (/ms) 
  b1p = -0.0316785972189014     (/mV) 
  c1p = 0.00043751852678292985     (/mV2) 
  d1p = -2.182210337922452e-06     (/mV3) 
  b2p = -0.025676454329437326     (/mV) 
  c2p = 0.0007017500668219226     (/mV2) 
  d2p = -3.860343090705927e-06     (/mV3) 

  am = 0.10607182980447505     (/mV) 
  bm = -3.320763613553543     (1) 
  vhm = -44.652400471232134     (mV) 
  Am = 0.15308898666019793     (/ms) 
  b1m = 0.09039234248432201     (/mV) 
  c1m = 0.0014357186171899483     (/mV2) 
  d1m = 1.727719641148694e-05     (/mV3) 
  b2m = 0.02218700036997514     (/mV) 
  c2m = 2.2192161922659285e-05     (/mV2) 
  d2m = -5.476921532631914e-07     (/mV3) 

  ac = 0.10607232812065745     (/mV) 
  bc = -1.5175473413162608     (1) 
  vhc = -20.913671345096958     (mV) 
  Ac = 0.3540896067295243     (/ms) 
  b1c = -0.04142225011605694     (/mV) 
  c1c = 0.0002820135971321023     (/mV2) 
  d1c = -8.164004575153624e-07     (/mV3) 
  b2c = -0.0757190852018145     (/mV) 
  c2c = -0.0005635664635673495     (/mV2) 
  d2c = -3.7987854749713668e-06     (/mV3) 

  aq = -0.14171763682482902     (/mV) 
  bq = 9.566447873142081     (1) 
  vhq = -72.30494496009798     (mV) 
  Aq = 29.796733263179966     (/ms) 
  b1q = -0.0918103704375076     (/mV) 
  c1q = 0.0008866713831835767     (/mV2) 
  d1q = -2.659763349106503e-06     (/mV3) 
  b2q = -0.09570774106920192     (/mV) 
  c2q = -0.0025256298395343667     (/mV2) 
  d2q = -4.608039612293257e-05     (/mV3) 
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
  hInf 
  hTau 
  pInf 
  pTau 
  mInf 
  mTau 
  cInf 
  cTau 
  qInf 
  qTau 
}

STATE
{
  n
  h
  p
  m
  c
  q
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*n*n*n*h*p*p*p*m*m*m*c*c*c*q
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
  h' = (hInf - h) / hTau 
  p' = (pInf - p) / pTau 
  m' = (mInf - m) / mTau 
  c' = (cInf - c) / cTau 
  q' = (qInf - q) / qTau 
}

INITIAL
{
  rates(v)
  n = nInf 
  h = hInf 
  p = pInf 
  m = mInf 
  c = cInf 
  q = qInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    cInf = 1/(1 + exp(-ac*v + bc)) 
    cTau = Ac / ( exp(-(b1c*(v-vhc) + c1c*(v-vhc)^2 + d1c*(v-vhc)^3)) + exp((b2c*(v-vhc) + c2c*(v-vhc)^2 + d2c*(v-vhc)^3)) ) 

    qInf = 1/(1 + exp(-aq*v + bq)) 
    qTau = Aq / ( exp(-(b1q*(v-vhq) + c1q*(v-vhq)^2 + d1q*(v-vhq)^3)) + exp((b2q*(v-vhq) + c2q*(v-vhq)^2 + d2q*(v-vhq)^3)) ) 


  UNITSON
}