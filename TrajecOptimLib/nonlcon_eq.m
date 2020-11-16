function ceq = nonlcon_eq(X,constants)

    X = matrify(X(:),constants);

    ceq = [];

    % dynamics
    ceq_dyn = cell(constants.Nv,1);
    for i = 1:constants.Nv
        [~,ceq_i] = constants.dynamics(X(:,:,i),constants);
        ceq_dyn{i} = ceq_i;
    end
    ceq = [ceq; cell2mat(ceq_dyn)];

end