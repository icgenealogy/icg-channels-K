NEURON {
    SUFFIX kdr
}
NEURON {
    USEION k READ ek WRITE ik
}
ASSIGNED {
    ik
    ek (mV)
}
PARAMETER {
	:erev 		= -90    (mV)
	gmax 		= 0.08   (umho)
        vrest           = 0

	mvalence 	= 3.
	mgamma 		= 0.7
	mbaserate 	= 0.1
	mvhalf 		= -35.
	mbasetau 	= 1.
	mtemp 		= 24.
        mq10            = 3
	mexp 		= 4

	hvalence 	= 0
	hgamma		= 0
	hbaserate 	= 0
	hvhalf 		= 0
	hbasetau 	= 0
	htemp 		= 0
        hq10            = 3
	hexp 		= 0

	cao (mM)
	cai (mM)
	celsius (degC)
	dt (ms)
	v (mV)

	vmax 		= 50   (mV)
	vmin 		= -100 (mV)
} : end PARAMETER

INCLUDE "custom_code/inc_files/64229_bg_cvode.inc"
PROCEDURE iassign() { i = g * (v - ek) ik=i }

:** kdrp
