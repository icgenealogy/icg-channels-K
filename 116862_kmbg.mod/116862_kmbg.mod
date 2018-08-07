NEURON {
    SUFFIX kmbg
}
NEURON {  USEION k READ ek WRITE ik }
ASSIGNED {
	ik
	ek (mV)
}
PARAMETER {
	:erev 		= -90    (mV)
	:ek		(mV)
        gbar 		= 0.1    (S/cm2)
        vrest           = 0.

	mvalence 	= 3.5
	mgamma 		= 0.8
	mbaserate 	= 0.2
	mvhalf 		= -46.
	mbasetau 	= 3.
	mtemp 		= 37
	mq10		= 3
	mexp 		= 3

	hvalence 	= 0
	hgamma		= 0
	hbaserate 	= 0
	hvhalf 		= 0
	hbasetau 	= 0
	htemp 		= 0
	hq10		= 3
	hexp 		= 0


	celsius			   (degC)
	dt 				   (ms)
	v 			       (mV)

	vmax 		= 100  (mV)
	vmin 		= -100 (mV)
} : end PARAMETER

INCLUDE "custom_code/inc_files/116862_bg_cvode.inc"
PROCEDURE iassign () { i = g*(v-ek) ik=i }
:* SYNAPSES
:** AMPA
