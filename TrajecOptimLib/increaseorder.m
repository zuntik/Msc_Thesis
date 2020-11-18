function Xout = increaseorder(X, constants, neworder)

    X = matrify(X, constants);

    X = BernsteinDegrElev(X, neworder);
    
    Xstatevars = X(2:neworder,1:constants.numvars,:);
    Xinputs = X(:,constants.numvars+1:end,:);
    
    Xout = zeros((neworder-1)*constants.numvars+(neworder+1)*constants.numinputs, constants.Nv);
    for i = 1:constants.Nv
        Xout(:,i) = [reshape(Xstatevars(:,:,i),[],1);reshape(Xinputs(:,:,i),[],1)];
    end
    Xout = Xout(:);
    
end