TITLE Ohmic K Current

: Implemented in Rubin and Cleland (2006) J Neurophysiology

NEURON {
	SUFFIX kO
	:USEION O READ eO WRITE iO VALENCE 1
	USEION k READ ek WRITE ik
        RANGE gbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gbar = 1.0 (mho/cm2)
	:eO = -90 (mV)
}

ASSIGNED {
	ik (mA/cm2)
        ek (mV)
}

BREAKPOINT {
	ik = gbar*(v - ek)
}



