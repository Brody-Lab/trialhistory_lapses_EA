function[history_params] = make_history_params(etaC, etaE, tauC, tauE)

%{
This function takes the four history_related params, and organises them
into a struct.
etaC: magnitude of bias following a correct rightward choice
etaE: magnitude of bias following an incorrect rightward choice
tauC: exponential decay param for correct trials (between 0 and 1)
tauE: exponential decay param for incorrect trials (between 0 and 1)

%}


history_params.etaC = etaC;
history_params.etaE = etaE;
history_params.tauC = tauC;
history_params.tauE = tauE;
end