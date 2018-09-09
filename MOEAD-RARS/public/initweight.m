function W=initweight(objDim, N, ifadjust)
% init_weights function initialize a pupulation of subproblems structure
% with the generated decomposition weight and the neighborhood relationship.
%initialize direction vectors with unity simplex;
%compute their neighborhood relationship and get their corresponding weight vectors
global params;
    %to load the parameter from the given file.
    filename = sprintf('weight/W%uD_%u.dat',objDim,params.popsize);
    path('../weight',path);
    if (exist(filename, 'file'))
    % copy from the exist file. 
     allws = importdata(filename);
     W = allws';
    end
     W((W < 1.0E-5))  = 1.0E-5;
    
     W=(1./W)./repmat(sum(1./W),[objDim 1]);
end
