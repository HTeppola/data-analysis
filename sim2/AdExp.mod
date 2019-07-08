: AdExp IF (Adaptive Exponential Intgrate and Fire)
: All printfs are for debugging purposes and can be ignored


NEURON {

    POINT_PROCESS AdExp
	RANGE num
    RANGE G_l, Delta_T, tau_w, a, b, E_l, I_ext, C
    RANGE taue, taui, erev, irev
    RANGE V_reset, V_thre, tRISE,tot
    RANGE gRISEinh,gRISEexc, iEXT,tauINH,tauEXC,tRISE,tRISEnmda,tDECAYnmda,gAMPA,gGABA
    RANGE iNMDA,iAMPA,iGABA,iTotal,gGABA,gAMPA,gNMDA,mAMPA,mNMDA,rpeso,Erev,v0_block,k_block,tau_r, tau_d, maxcurrent
    RANGE label

}

UNITS {

	(mV) = (millivolt)
	(pA) = (picoamp)
	(uS) = (microsiemens)
    (nS) = (nanosiemens)
	(pS) = (picosiemens)
}

PARAMETER {
	num = -1  					:Cell ID
	V_reset = -60	(mV)		:reset potential after spike
	V_thre = -50.4	(mV)		:threshold for spike detection
	V_spike = -10	(mV)	    :value of spike
	a = 4			(nS)		:coupling with adaptive variable
	b = 80.5		(pA)		:adaptive increment
	tau_w = 144		(ms)		:adaptive time costant
	E_l = -70.6		(mV)		:resting potential for leak term
	G_l = 30		(nS)	    :leak conduptance
	Delta_T = 2		(mV)	    :speed of exp

	eps = 0.1		(ms) 		:smaller than time step
	iEXT = 0      (pA) 	    :external current
	C = 281 		(nS)		:membrane conduptance
    tauEXC = 3 		(ms) 		:excitatory time costant
    tauINH = 8		(ms) 		:inhibitory time costant
    erev = 0        (mV) 		:reverse potential for exc
    irev = -70		(mV)		:reverse potential for inh
    tRISE = 1 		(ms) 		:time to activate the synapse

    mNMDA=.05
    mAMPA=1
    v0_block 	= -50
	k_block 	= 8

	tRISEnmda = 5.63    (ms)        :Chapman DE 2003, Table 1 - rise tau
	tDECAYnmda = 140     (ms)        :Chapman DE 2003, Fig 2B Ventromedial - decay tau

	Erev = 0         (mV)        :reversal potential, Dalby 2003
	rpeso = 12        (nS)
    maxcurrent = 500 : maximum allowed total current
    label = 0 : assign an arbitrary label to the cell

}

ASSIGNED {
    gNMDA           (nS)
    gAMPA           (nS)
    gGABA           (nS)
	iAMPA           (pA)
    iGABA           (pA)
	iNMDA           (pA)
    iNOISE          (pA)
    iTotal          (pA)


}

STATE {
    vv              (mV)
    ww              (pA)
    gEXC            (nS)
    gINH            (nS)
    gRISEexc        (nS)
    gRISEinh        (nS)
    y1              (nS)
    y2              (nS)
}

INITIAL { :initializion and activation of the process

    gEXC=0          (nS)
    gINH=0          (nS)
    gRISEexc=0      (nS)
    gRISEinh=0      (nS)
    gAMPA=0         (nS)
    gGABA=0         (nS)
    ww=0            (pA)
    vv=-60          (mV)
    net_send(0,1)
}


BREAKPOINT {



SOLVE states METHOD derivimplicit


}

DERIVATIVE states {

	y1' = -y1 / (tDECAYnmda)
	y2' = -y2 / (tRISEnmda)
    gNMDA = (y2-y1)
	iNMDA = (gNMDA)* (vv - Erev)*1 / ( 1 + exp ( - ( vv - v0_block ) / k_block ))

    gRISEexc' = - gRISEexc/tRISE        :rise of excitatory synapse
    gRISEinh'= - gRISEinh/tRISE         :rise of inhibitory synapse
    gEXC' = - gEXC /tauEXC              :decay of exc conduptance
    gINH' = - gINH /tauINH              :decay of inh conduptance

    gAMPA = (gRISEexc-gEXC)
	iAMPA = gAMPA * ( vv - erev )

    gGABA = (gRISEinh-gINH)
    iGABA = -gGABA * ( vv - irev )


    														:iNOISE = - (x2-x1) * ( vv - erev )
    iTotal=-ww+iEXT+iAMPA+iGABA+iNMDA  						:+iNOISE
    if (iTotal>maxcurrent) {
        iTotal = maxcurrent
    }
    
    :if (vv > V_thre) {
    :  V_thre = vv+1.0
    :}

    vv'=1/C*(-G_l*(vv-E_l)+G_l*Delta_T*exp(((vv-V_thre))/Delta_T)+iTotal)
    ww' = 1/tau_w*(a*(vv-E_l)-ww)

    }


NET_RECEIVE (w0,w1,w2) {
    INITIAL {
        w1=w1
        w2=w2
        }
    if (flag == 1) {                            :wait for membrane potential to reach threshold
        WATCH (vv>V_spike) 2
    }
    else if (flag == 2) {                       :threshold reached, fire and reset the model

    net_event(t)
    vv = V_reset
    ww = ww+b

    }
    else if (flag==0) {                         :an event has arrived -> synaptic activation

    if ( w0 > 0 ) {                             :positive weight -> excitatory synapse

        if (w0>1000) {                              :overlimit conductance -> poisson process
            gEXC=gEXC+rpeso
            gRISEexc=gRISEexc+rpeso
            :x1=x1+rpeso
            :x2=x2+rpeso
        }

        else {                                      :genuine excitation -> update synconductance

            y1=y1+(w0 *(mNMDA))
            y2=y2+(w0 *(mNMDA))
            gEXC=gEXC+w0 *w1*(mAMPA)
            gRISEexc=gRISEexc+w0 *w1*(mAMPA)
        }

    }
    else {                                          :genuine inhibition -> update synconductance

        gINH=gINH+w0
        gRISEinh=gRISEinh+w0

        }
    }
}
