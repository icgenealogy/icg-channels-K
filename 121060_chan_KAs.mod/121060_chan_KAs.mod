TITLE KAs channel

NEURON {
	SUFFIX KAs
	USEION  k READ ek WRITE ik
	RANGE gmax, g, ik
	GLOBAL minf, mtau, hinf, htau
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	gmax = 9.51e-4	(mho/cm2)	<0,1e9> : for middle and distal dendrites
	ek	        (mV)
	m_vh = -27	(mV)	: half activation
	m_ve = -16	(mV)	: slope
	h_vh = -33.5	(mV) : half activation
	h_ve = 21.5	(mV) : slope
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	ik	(mA/cm2)
	minf	(1)
	hinf   (1)
	mtau	(ms)
	htau  (ms)
}

CONSTANT {
	a = 0.996 (1)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*m*m*(a*h+(1-a))
	ik = g*(v - ek)
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states { 
	rates(v)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
}

: rates() computes rate and other constants at present v
: call once from hoc to initialize inf at resting v
PROCEDURE rates(v(mV)) {
	LOCAL m_alpha, m_beta
	TABLE mtau, minf, htau, hinf FROM -120 TO 40 WITH 160
UNITSOFF
	: "m" KAs activation and "h" KAs inactivation
	mtau = 0.378+9.91*exp(-((v+34.3)/30.1)^2)
	minf = 1/(1 + exp((v - m_vh)/m_ve))
	
	m_alpha = exp(-(v+90.96)/29.01)
	m_beta = exp((v+90.96)/100)
	htau = 1097.4/(m_alpha+m_beta)
	hinf = 1/(1 + exp((v - h_vh)/h_ve))	
}
UNITSON
