function [Aeq, beq] = LinearConstr_fd(n,T,pos0,vel0)

    % impose p' = v <=> -p'+v=0

    Aeq_dyn = [ BernsteinDerivMat(n+1,T)*BernsteinDegrElevMat(n,n+1), -eye(n+1), zeros(n+1)];
    beq_dyn = zeros(n+1,1);

    % impose initial conditions and final conditions

    Aeq_initial = zeros(3,3*n+3);
    Aeq_initial(1,1) = 1;
    Aeq_initial(2,1) = -n/T;
    Aeq_initial(2,2) = n/T;
    Aeq_initial(3,n+1)=1;
    beq_initial = [ pos0 ; vel0; 0];

    Aeq = [ Aeq_dyn ; Aeq_initial ];
    beq = [ beq_dyn ; beq_initial ];

end