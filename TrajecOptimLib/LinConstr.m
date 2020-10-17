function [Aeq, beq] = LinConstr(constants)

    N = constants.N;
    numvars = constants.numvars;
    numinputs = constants.numinputs;
    Nv = constants.Nv;
    xi = constants.xi;
    xf = constants.xf;

    Aeq = zeros(Nv*numvars*2, (N+1)*(numvars+numinputs)*Nv);
    beq = zeros(Nv*numinputs*2,1);
    for veh = 1:Nv
        for inp = 1:numvars
            a=zeros(N+1, numvars+numinputs, Nv);
            b=zeros(N+1, numvars+numinputs, Nv);
            a(1, inp, veh) = 1;
            b(N+1, inp, veh) = 1;
            Aeq((veh*(numvars-1)+inp)*2-1,:)=a(:);
            Aeq((veh*(numvars-1)+inp)*2,:)=b(:);
            beq((veh*(numvars-1)+inp)*2-1)=xi(veh,inp);
            beq((veh*(numvars-1)+inp)*2)=xf(veh,inp);
        end
    end

end
