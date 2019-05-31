NEURON {
	SUFFIX ionleak
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gk, gna, ik, ina
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	gk = 0 (S/cm2)
	gna = 0 (S/cm2)
}

ASSIGNED {
	v (mV)
	ena (mV)
	ek (mV)
	ik (mA/cm2)
	ina (mA/cm2)
}

INITIAL {
}

BREAKPOINT {
	ina = gna*(v - ena)
	ik = gk*(v - ek)
}

