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
	:erev 		= -75.  (mV)
	gbar 		= 0.015    (mho/cm2)


        vrest           = -60
	exptemp		= 37
	maflag 		= 3
	malphaA 	= -0.016
	malphaB		= -5.0
	malphaV0	= 35.1
	mbflag 		= 1
	mbetaA 		= 0.25
	mbetaB		= -40.
	mbetaV0		= 20.
	mq10		= 3
	mexp 		= 1

	haflag 		= 0
	halphaA 	= 0
	halphaB		= 0
	halphaV0	= 0
	hbflag 		= 0
	hbetaA 		= 0
	hbetaB		= 0
	hbetaV0		= 0
	hq10		= 3
	hexp 		= 0

	vmax 		= 100  (mV)
	vmin 		= -100 (mV)
} : end PARAMETER

INCLUDE "custom_code/inc_files/35358_geneval_cvode.inc"

PROCEDURE iassign () { i = g*(v-ek) ik=i }
:* cal from calRT03 in ~/nrniv/place/mod/parameters.multi
