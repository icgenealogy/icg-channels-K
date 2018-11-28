TITLE KRP channel
COMMENT
Unit check passed even without Units off
ENDCOMMENT

NEURON {
	SUFFIX KRP
	USEION  k READ ek WRITE ik
	RANGE g, gmax, ik
	GLOBAL minf, mtau, hinf, htau
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	gmax = 0.001	(mho/cm2)	<0,1e9>
	ek 		        (mV)
	m_vh = -13.5	(mV)	: half activation
	m_ve = -11.8	(mV)	: slope
	h_vh = -54.7	(mV) : half activation
	h_ve = 18.6	(mV) : slope
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
	a = 0.7 (1)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*m*(a*h+(1-a))
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

FUNCTION_TABLE tabmtau(v(mV)) (ms)
FUNCTION_TABLE tabhtau(v(mV)) (ms)

: rates() computes rate and other constants at present v
: call once from hoc to initialize inf at resting v
PROCEDURE rates(v(mV)) {
:	TABLE minf, hinf FROM -120 TO 30 WITH 70
	mtau = tabmtau(v)
	minf = 1/(1 + exp((v - m_vh)/m_ve))
	htau = tabhtau(v)
	hinf = 1/(1 + exp((v - h_vh)/h_ve))
}

