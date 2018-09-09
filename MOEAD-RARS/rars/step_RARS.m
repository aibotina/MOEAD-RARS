function step_RARS(mop)

global params population data_r_nei_ sum_iter maxfes;

    subindex = prob_select();
    nsub     = length(subindex);
    subindex = subindex(randperm(nsub));
    for isub = subindex
        % increase FES of each sub-problem
        params.count(isub) = params.count(isub) + 1;
        % dynamic neighborhood
        if rand() < params.pns
            neighborhood = population.neighbor(1:params.niche,isub);
        else
            neighborhood = 1:params.popsize;
        end

        % new point generation using genetic operations, and evaluate it.
        ind     = de(isub, neighborhood);
        obj     = evaluate(mop, ind);
        
        % update the idealpoint.
        population.ideapoint  = min(population.ideapoint, obj);

        % update the neighbours.
        update(obj, ind, 1:params.popsize,isub);
        sum_iter(isub)=sum_iter(isub)+1;
    end
    
    % increase FES
    params.fes  = params.fes + nsub;
    % update iterate counter
    params.gen  = params.gen + 1;
    
     utility_update();
    % save current population
    if params.S(1) == '2' || params.S(3) == '2' || params.S(3) == '3'
        population.hp(1:(2*params.dt-1)*params.fdim,:)       = population.hp((params.fdim+1):end,:);
        population.hp((2*params.dt-1)*params.fdim+1:end,:)   = population.objective(:,:);
         for i=1:20
             data_r_nei_{i}=data_r_nei_{i+1};
         end
    end
    clear subindex rndi ind obj neighbourhood;
end



%% update strategy
function update(obj, parameter, neighborhood,isubi)
%updating of the neighborindex with the given new individuals
global params population data_r_nei_ sum_flag_iter maxfes;

    % objective values of current solution
    newobj  = subobjective(population.W(:,neighborhood), obj, population.ideapoint, params.dmethod);
    % previous objective values
    oldobj  = subobjective(population.W(:,neighborhood), population.objective(:,neighborhood), population.ideapoint, params.dmethod);    
    % new solution is better?
    improve = (oldobj - newobj) ./ oldobj;
    [iv,id] = max(improve);
       if iv > 0 
         population.parameter(:,id) = parameter;
         population.objective(:,id) = obj;  
         population.aindex(id)      = 1;
         data_r_nei_{21}(isubi,id)= data_r_nei_{21}(isubi,id)+iv;
         sum_flag_iter(isubi)=sum_flag_iter(isubi)+1;
       end
   
    clear subp newobj oops oldobj C toupdate;
end

%%
function ind = de(index, neighborhood)

global population params maxfes;
    %parents
    si      = ones(1,3)*index;
    c       = 0;
    while si(2)==si(1) || si(3)==si(1) || si(3)==si(2)
        si(2)   = neighborhood(floor(rand()*length(neighborhood))+1);
        si(3)   = neighborhood(floor(rand()*length(neighborhood))+1);
        c       = c + 1;
        if c>10, disp('error: de'); break; end % to ensure there is no dead loop
    end
  selectpoints    = population.parameter(:, si);

    %generate new trial point  CR ²ÎÊý¿ØÖÆ
    %  UF MOP
   newpoint        = selectpoints(:,1)+params.F*(selectpoints(:,2)-selectpoints(:,3));
   %WFG
%    if rand()<0.5 
%        newpoint        = selectpoints(:,1)+params.F*(selectpoints(:,2)-selectpoints(:,3));
%    else
%        newpoint=selectpoints(:,1);
%    end

    %repair the new value
    rnds            = rand(params.xdim,1);
    pos             = newpoint>params.xupp;
    if sum(pos)>0
        newpoint(pos) = selectpoints(pos,1) + rnds(pos,1).*(params.xupp(pos)-selectpoints(pos,1));
    end
    pos             = newpoint<params.xlow;
    if sum(pos)>0
        newpoint(pos) = selectpoints(pos,1) - rnds(pos,1).*(selectpoints(pos,1)-params.xlow(pos));
    end
    
    ind             = polymutate(newpoint, params.pm);

    clear si selectpoints newpoint pos;
end

%%
function ind = polymutate(ind, rate)
global params;

    eta_m   = params.etam;
    mut_pow = 1.0 / (eta_m + 1.0);   
    for j = 1:params.xdim
      r = rand();
      if r <= rate
        y       = ind(j);
        yl      = params.xlow(j);
        yu      = params.xupp(j);
        delta1  = (y - yl) / (yu - yl);
        delta2  = (yu - y) / (yu - yl);
        
        rnd     = rand();
        if (rnd <= 0.5) 
          xy    = 1.0 - delta1;
          val   = 2.0 * rnd + (1.0 - 2.0 * rnd) * (xy^(eta_m + 1.0));
          deltaq= (val^mut_pow) - 1.0;
        else 
          xy    = 1.0 - delta2;
          val   = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (xy^ (eta_m + 1.0));
          deltaq= 1.0 - (val^mut_pow);
        end

        y   = y + deltaq * (yu - yl);        

%         if y < yl, y = yl; end
%         if y > yu, y = yu; end
%         ind(j) = y;        
        if y >= yl && y <= yu
            ind(j) = y;
        end
      end
    end
end

%% select sub-problems to evolve
function index = prob_select()

global  params
 
if params.S(1) == '0' 
	index   = 1:params.popsize;
else
    ind     = rand(1,params.popsize) <= params.ps;
    index   = find(ind > 0);
end

% it is not likey to happen, just in case
if isempty(index)
    index   = 1:params.popsize;
end

clear ind;
end

%% utility function update
function utility_update()

global population params;

% step 1: utility function
switch params.S(1)
    case {'0','1'}  % no update
        return;
    case '2'        % improvement
        v = utility_improvement;                
end

% record the utility value
params.hvs(1:(params.dt-1),:)    = params.hvs(2:end,:);
params.hvs(end,:)                = v;

% step 2: utility function aggregation
switch params.S(2)
    case '0'    % no aggreation
    case '1'    % neighborhood aggregation
        v = zeros(1, params.popsize);
        for i=1:params.popsize
            nb      = population.neighbor(1:params.niche,i);
            vt      = params.hvs(end,nb);
            v(i)    = mean(vt(:));
        end
        v = v / max(v);
        clear nb vt;        
    case '2'    % history aggregation
        v = mean(params.hvs);
        v = v / max(v);
    case '3'    % both neighborhood and history aggregation
        v = zeros(1, params.popsize);
        for i=1:params.popsize
            nb      = population.neighbor(1:params.niche,i);
            vt      = params.hvs(:,nb);
            v(i)    = mean(vt(:));
        end
        v = v / max(v);
        clear nb vt;        
end

% output
params.ps = v;
clear v;

end
%%
function v = utility_improvement()
    global  data_r_nei_ cc population maxfes params iternum;
    r_nei_= data_r_nei_{21}-data_r_nei_{1};
    d_r_nei=diag(r_nei_);
    aaa=0.8;
    sum_row=sum(r_nei_);
    sum_cow=sum(r_nei_,2)';
    for i=1:size(r_nei_,1)  
         v(i)=(aaa*sum_row(i)+(1-aaa)*sum_cow(i)+1.0E-50)/(max(aaa*sum_row+(1-aaa)*sum_cow)+1.0E-50);   
    end
    if  params.fes<(0.05)*maxfes
    %   nei size
    ns_gg=15;
     if size(population.W,1)==2
             [vv1,ss1]=sort(population.W(1,:));
             [vv2,ss2]=sort(population.W(2,:));
             nei_gg1 = population.neighbor(1:params.niche,ss1(1));
            nei_gg2 = population.neighbor(1:params.niche,ss2(1));
            flag_n1=0;
            flag_n2=0;
            for i=1:ns_gg
                if v(nei_gg1(i))>0.5
                    flag_n1=1;
                    break;
                end
            end
            for i=1:ns_gg
                if v(nei_gg2(i))>0.5
                    flag_n2=1;
                    break;
                end
            end
            if flag_n1==0
                v(nei_gg1(floor(rand()*ns_gg)+1))=1;
            end
            if flag_n2==0
                v(nei_gg2(floor(rand()*ns_gg)+1))=1;
            end
            
         else
            [vv1,ss1]=sort(population.W(1,:));
            [vv2,ss2]=sort(population.W(2,:));
            [vv3,ss3]=sort(population.W(3,:));
            nei_gg1 = population.neighbor(1:params.niche,ss1(1));
            nei_gg2 = population.neighbor(1:params.niche,ss2(1));
            nei_gg3 = population.neighbor(1:params.niche,ss3(1));  
            flag_n1=0;
            flag_n2=0;
            flag_n3=0;
            for i=1:ns_gg
                if v(nei_gg1(i))>0.5
                    flag_n1=1;
                    break;
                end
            end
            for i=1:ns_gg
                if v(nei_gg2(i))>0.5
                    flag_n2=1;
                    break;
                end
            end
             for i=1:ns_gg
                if v(nei_gg3(i))>0.5
                    flag_n3=1;
                    break;
                end
            end
            if flag_n1==0
                v(nei_gg1(floor(rand()*ns_gg)+1))=1;
            end
            if flag_n2==0
                v(nei_gg2(floor(rand()*ns_gg)+1))=1;
            end
             if flag_n3==0
                v(nei_gg3(floor(rand()*ns_gg)+1))=1;
             end
     end
    end
    clear cpop hpop cobj hobj;    
end








