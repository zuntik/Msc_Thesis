function xinit = rand_init_guess(constants)
    xinit = rand((constants.numvars*(constants.N-1)+...
        constants.numinputs*(constants.N+1))*constants.Nv,1);
end