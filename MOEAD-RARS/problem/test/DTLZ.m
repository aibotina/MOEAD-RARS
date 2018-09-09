function mop=DTLZ(testname,pdim,odim)
%%run for the test problems dtlz serious.
mop=struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);
switch testname
    case 'DTLZ1'
        mop=dtlz1(mop,pdim,odim);
    case 'DTLZ2'
        mop=dtlz2(mop,pdim,odim);        
    case 'DTLZ3'
        mop=dtlz3(mop,pdim,odim);        
    case 'DTLZ4'
        mop=dtlz4(mop,pdim,odim);        
    case 'DTLZ5'
        mop=dtlz5(mop,pdim,odim);        
    case 'DTLZ6'
        mop=dtlz6(mop,pdim,odim);        
    case 'DTLZ7'
        mop=dtlz7(mop,pdim,odim);                                                
    otherwise
        error('Undefined test problem name');
end
%%%%%%%%%%FUNCTIONS%%%%%%
    function p=dtlz1(p,pdim,odim)
        p.name='DTLZ1';
        p.od=odim;
        p.pd=pdim;
        p.domain=[zeros(pdim,1) ones(pdim,1)];
        p.func=@evaluate;
        function y=evaluate(x)
            k=pdim-odim+1;
            y=zeros(odim,1);
            xm=x(pdim-k+1:pdim,1);
            g=100*(k+sum((xm-0.5).^2-cos(20*pi*(xm-0.5))));
            y(1)=0.5*prod(x(1:odim-1,1))*(1+g);
            y(odim)=0.5*(1-x(1))*(1+g);
            for i=2:odim-1
                y(i)=0.5*prod(x(1:odim-i,1))*(1-x(odim-i+1))*(1+g);
            end                        
        end
    end
%%%%%%%
    function p=dtlz2(p,pdim,odim)
        p.name='DTLZ2';
        p.od=odim;
        p.pd=pdim;
        p.domain=[zeros(pdim,1) ones(pdim,1)];
        p.func=@evaluate;
        function y=evaluate(x)
            k=pdim-odim+1;
            y=zeros(odim,1);
            xm=x(pdim-k+1:pdim,1);
            g=sum((xm-0.5).^2);
            y(1)=prod(cos(pi*x(1:odim-1,1)./2))*(1+g);
            y(odim)=sin(pi*x(1)/2)*(1+g);
            for i=2:odim-1
                y(i)=prod(cos(pi*x(1:odim-i,1)./2))*sin(pi*x(odim-i+1)./2)*(1+g);
            end  
        end
    end
%%%%%%%
    function p=dtlz3(p,pdim,odim)
        p.name='DTLZ3';
        p.od=odim;
        p.pd=pdim;
        p.domain=[zeros(pdim,1) ones(pdim,1)];
        p.func=@evaluate;
        function y=evaluate(x)
            k=pdim-odim+1;
            y=zeros(odim,1);
            xm=x(pdim-k+1:pdim,1);
            g=100*(k+sum((xm-0.5).^2-cos(20*pi*(xm-0.5))));
            y(1)=prod(cos(pi*x(1:odim-1,1)./2))*(1+g);
            y(odim)=sin(pi*x(1)/2)*(1+g);
            for i=2:odim-1
                y(i)=prod(cos(pi*x(1:odim-i,1)./2))*sin(pi*x(odim-i+1)./2)*(1+g);
            end  
        end
    end
%%%%%%%
    function p=dtlz4(p,pdim,odim)
        p.name='DTLZ4';
        p.od=odim;
        p.pd=pdim;
        p.domain=[zeros(pdim,1) ones(pdim,1)];
        p.func=@evaluate;
        function y=evaluate(x)
            a=100;
            x=x.^a;
            k=pdim-odim+1;
            y=zeros(odim,1);
            xm=x(pdim-k+1:pdim,1);
            g=sum((xm-0.5).^2);
            y(1)=prod(cos(pi*x(1:odim-1,1)./2))*(1+g);
            y(odim)=sin(pi*x(1)/2)*(1+g);
            for i=2:odim-1
                y(i)=prod(cos(pi*x(1:odim-i,1)./2))*sin(pi*x(odim-i+1)./2)*(1+g);
            end  
        end
    end
%%%%%%%
    function p=dtlz5(p,pdim,odim)
        p.name='DTLZ5';
        p.od=odim;
        p.pd=pdim;
        p.domain=[zeros(pdim,1) ones(pdim,1)];
        p.func=@evaluate;
        function y=evaluate(x)
            k=pdim-odim+1;
            y=zeros(odim,1);
            xm=x(pdim-k+1:pdim,1);
            g=sum(xm.^0.1);
            y(1)=prod(cos(pi*x(1:odim-1,1)./2))*(1+g);
            y(odim)=sin(pi*x(1)/2)*(1+g);
            zta=(pi/(4*(1+g)))*(1+2*g*x);
            for i=2:odim-1
                y(i)=prod(cos(zta(i)))*sin(zta(i+1))*(1+g);
            end  
        end
    end
%%%%%%%
    function p=dtlz6(p,pdim,odim)
        p.name='DTLZ6';
        p.od=odim;
        p.pd=pdim;
        p.domain=[zeros(pdim,1) ones(pdim,1)];
        p.func=@evaluate;
        function y=evaluate(x)
            k=pdim-odim+1;
            y=zeros(odim,1);
            xm=x(pdim-k+1:pdim,1);
            for i=1:odim-1
                y(i)=x(i);
            end
            g=1+(9/k)*sum(xm);
            h=odim-sum((y(1:odim-1)./(1+g)).*(1+sin(3*pi*y(1:odim-1))));
            y(odim)=(1+g).*h;
            
        end
    end
%%%%%%%
    function p=dtlz7(p,pdim,odim)
        p.name='DTLZ7';
        p.od=odim;
        p.pd=pdim;
        p.domain=[zeros(pdim,1) ones(pdim,1)];
        p.func=@evaluate;
        function y=evaluate(x)
            y=zeros(odim,1);
            g=zeros(odim,1);
            for i=1:odim
                y(i)=(1/(floor(pdim/odim)))*sum(x((floor((i-1)*(pdim/odim))+1):(floor(i*pdim/odim))),1);
            end
            for j=1:odim-1
                g(j)=y(odim)+4*y(j)-1;
            end
            k=1;
            for i=1:odim-1
                for j=1:odim-1
                    if i~=j
                        z(k)=y(i)+y(j);
                        k=k+1;
                    end
                end
            end
            g(odim)=2*y(odim)+min(z)-1;
            k1=1;
            for i=1:odim
                if g(i)<0
                    k1=k1+1;
                end
            end
            if k1>=2
                y=y+1000;
            end
        end
    end
%%%%%%%






end

