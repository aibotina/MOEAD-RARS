function demo_RARS()

global params maxfes population data_r_nei_ PFStar sum_flag_iter sum_iter flag_r_nei_ iternum;

path('../problem',path); 
path('../problem/cec09',path); 
path('../problem/test',path); 
path('./PFStar',path); 
path('../public',path);
re_pf_RARS=cell(9,30);
data_r_nei_=cell(1,21);
name_func={'tec09_f2'};
for seq=1:1
    for nrun=1:2
    problem   = char(name_func(seq));
    iternum=0;
   mop     = testmop(problem,30);
   popsize =300;
   maxfes  = 150000; 
PFStar   = load(strcat('PFStar/',problem,'.dat'));
tic;
init('problem', mop, 'popsize', popsize, 'niche', 20, 'pns', 0.8, 'F', 0.5, 'S', '20', 'DT', 20, 'method','ts');
for i=1:21
  data_r_nei_{i}=zeros(params.popsize);
end
  sum_iter=zeros(params.popsize,1);
  sum_flag_iter=zeros(params.popsize,1);
  flag_r_nei_=zeros(params.popsize);
g = 1;
while params.fes < maxfes
    step_RARS(mop);
    g =  g+1;
end
endt    = toc;

PF=population.objective;
re_pf_RARS{seq,nrun}=PF;
disp(endt); 
    end
end
save re_pf_RARS;
end
