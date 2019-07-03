COMMENT

A "simple" implementation of the Izhikevich neuron with AMPA, NMDA,
GABA_A, and GABA_B receptor dynamics. Parameter values are taken
from the appendix to
   Izhikevich EM, Edelman GM (2008).
   "Large-scale model of mammalian thalamocortical systems." 
   PNAS 105(9) 3593-3598.

Example usage (in Python):
  from neuron import h
  dummycell = h.Section() # Since Izhi is a point process, it needs to be in a section
  cells = [0,1] # Initialize cells
  for c in cells: cells[c] = h.Izhi(0) # Create two new Izhikevich cells
  connection = h.NetCon(cells[0], cells[1]) # Connect them
  
See also (if available) izhi.py, a Python wrapper for different types of
Izhikevich neuron.

Version: 2013oct16 by Cliff Kerr (cliffk@neurosim.downstate.edu)

ENDCOMMENT


: Declare name of object and variables
NEURON {
  POINT_PROCESS Izhi
  RANGE C, k, vr, vt, vpeak, a, b, c, d, tauAMPA, tauNMDA, tauGABAA, tauGABAB, celltype, alive, cellid, verbose
  RANGE V, u, t0, I
}

: Specify units that have physiological interpretations (NB: ms is already declared)
UNITS {
	(mV) = (millivolt)
	(uM) = (micrometer)
}


: Parameters from Izhikevich 2008, PNAS for P2/3 cell, no dendrites
PARAMETER {
  C = 100 : Capacitance
  k = 3
  vr = -60 (mV) : Resting membrane potential
  vt = -50 (mV) : Membrane threhsold
  vpeak = 30 (mV) : Peak voltage
  a = 0.01
  b = 5
  c = -60
  d = 400
  tauAMPA = 5 (ms) : Receptor time constant, AMPA
  tauNMDA = 150 (ms) : Receptor time constant, NMDA
  tauGABAA = 6 (ms) : Receptor time constant, GABAA
  tauGABAB = 150 (ms) : Receptor time constant, GABAB
  celltype = 1 : A flag for indicating what kind of cell it is,  used for changing the dynamics slightly. Cell types are: 1 = pyramidal, 2 = fast spiking, 3 = low threshold, 4 = thalamocortical, 5 = reticular
  alive = 1 : A flag for deciding whether or not the cell is alive -- if it's dead, acts normally except it doesn't fire spikes
  cellid = -1 : A parameter for storing the cell ID, if required (useful for diagnostic information)
  verbose = 0 : Whether or not to print diagnostic information to file -- WARNING, do not modify this manually -- it's set by useverbose()
}


: Variables used for internal calculations
ASSIGNED {
  t0 : Previous time
  factor : Voltage factor used for calculating the current
  eventflag : For diagnostic information
  V (mV) : Membrane voltage
  u (mV) : Slow current/recovery variable
  gAMPA : AMPA conductance
  gNMDA : NMDA conductance
  gGABAA : GABAA conductance
  gGABAB : GABAB conductance
  I : Total current
  delta : Time step
}


: Initial conditions -- from, like, every Izhikevich paper
INITIAL {
  V = vr
  u = 0.2*V
  t0 = t
  gAMPA = 0
  gNMDA = 0
  gGABAA = 0
  gGABAB = 0
  I = 0
  delta = 0
  net_send(0,1) : Required for the WATCH statement to be active
}


: Function for printing diagnostic information to a file -- usage example: cell.useverbose(2,"logfile.txt")
VERBATIM
char filename[1000]; // Allocate some memory for the filename
ENDVERBATIM
PROCEDURE useverbose() { : Create user-accessible function
  VERBATIM
  #include<stdio.h> // Basic input-output
  verbose = (float) *getarg(1); // Set verbosity -- 0 = none, 1 = events, 2 = events + timesteps
  strcpy(filename, gargstr(2)); // Copy input filename into memory
  ENDVERBATIM
}



: Define neuron dynamics
BREAKPOINT {
  delta = t-t0 : Find time difference

  : Receptor dynamics -- the correct form is gAMPA = gAMPA*exp(-delta/tauAMPA), but this is 30% slower and, in the end, not really any more physiologically realistic
  gAMPA = gAMPA - delta*gAMPA/tauAMPA : "Exponential" decays -- fast excitatory (AMPA)
  gNMDA = gNMDA - delta*gNMDA/tauNMDA : Slow excitatory (NMDA)
  gGABAA = gGABAA - delta*gGABAA/tauGABAA : Fast inhibitory (GABA_A)
  gGABAB = gGABAB - delta*gGABAB/tauGABAB : Slow inhibitory (GABA_B)
  
  : Calculate current
  factor = ((V+80)/60)*((V+80)/60)
  I = gAMPA*(V-0) + gNMDA*factor/(1+factor)*(V-0) + gGABAA*(V+70) + gGABAB*(V+90)
  
  : Calculate neuronal dynamics; -I since I = -I_{syn}, which is really what I is as I've defined it above
  V = V + delta*(V + k*(V-vr)*(V-vt) - u - I)/C  : Calculate voltage
  u = u + delta*a*(b*(V-vr)-u) : Calculate recovery variable
  
  : Cell-type specific dynamics
  if (celltype>2) {
     : For LTS neurons, reset u if it gets too big
     if (celltype==3) {
       if (u>670) {u=670}
     }
     
     : For TCR neurons, reset b
     if (celltype==4) {
       if (V>-65) {b=0}
       else {b=15}
     }
     
     : For TRN neurons, reset b
     if (celltype==5) {
       if (V>-65) {b=2}
       else {b=10}
     }
  }

  t0=t : Reset last time so delta can be calculated in the next time step
  
  : Print diagnostic inormation to a file
  if (verbose>1) { : Verbose turned on?
    VERBATIM
    FILE *outfile; // Declare file object
    outfile=fopen(filename,"a"); // Open file for appending
    fprintf(outfile,"%8.2f   cell=%6.0f   delta=%8.2f   gAMPA=%8.2f   gNMDA=%8.2f   gGABAA=%8.2f   gGABAB=%8.2f   factor=%8.2f   I=%8.2f   V=%8.2f   u=%8.2f (timestep)\n",t,cellid,delta,gAMPA,gNMDA,gGABAA,gGABAB,factor,I,V,u);
    fclose(outfile); // Close file
    ENDVERBATIM
  }
  
  
}



: Input received
NET_RECEIVE (wAMPA, wNMDA, wGABAA, wGABAB) {  
  INITIAL { wAMPA=wAMPA wNMDA=wNMDA wGABAA=wGABAA wGABAB=wGABAB } : Insanely stupid but required, otherwise reset to 0, arrrgggh

  : Check if spike occurred
  if (flag == 1) { : Fake event from INITIAL block
    WATCH (V>vpeak) 2 : Check if threshold has been crossed, and if so, set flag=2
  } 
  
  : Event created by WATCH statement -- i.e. threshold crossed
  else if (flag == 2) { 
    if (alive) {net_event(t)} : Send spike event if the cell is alive
    V = c : Reset voltage
    u = u+d : Reset recovery variable
    gAMPA = 0 : Reset conductances -- not mentioned in Izhikevich's paper but necessary to stop things from exploding!
    gNMDA = 0
    gGABAA = 0
    gGABAB = 0
  } 
  
  : Actual input, calculate receptor dynamics
  else {
    gAMPA = gAMPA + wAMPA
    gNMDA = gNMDA + wNMDA
    gGABAA = gGABAA + wGABAA
    gGABAB = gGABAB + wGABAB
  }
  
  : Print diagnostic information to a file
  if (verbose>0) { : Verbose turned on?
    eventflag = flag
    VERBATIM
    FILE *outfile; // Declare file object
    outfile=fopen(filename,"a"); // Open file for appending
    fprintf(outfile,"%8.2f   cell=%6.0f   flag=%1.0f   gAMPA=%8.2f   gNMDA=%8.2f   gGABAA=%8.2f   gGABAB=%8.2f   V=%8.2f   u=%8.2f (event)\n",cellid,t,eventflag,gAMPA,gNMDA,gGABAA,gGABAB,V,u);
    fclose(outfile); // Close file
    ENDVERBATIM
  }
  
  
}
