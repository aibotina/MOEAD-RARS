function p=wfg1(p,odim,k,l)
        p.name='WFG1';
        p.od=odim;
        p.pd=k+l;
        if odim==2
            p.pf=load('WFG1.2D.pf');
        else
            p.pf=load('WFG1.3D.pf');
        end
        p.domain=[zeros(p.pd,1) 2.*(1:p.pd)'];
        p.func=@evaluate;
        function f=evaluate(z)
          %  z=z';
          
            cc=p.domain(:,:)';
           % z01=z./p.domain(2,:)';
           z01=z./cc(2,:)';
            y=z01;
            t1 = zeros(p.pd,1);t2 = zeros(p.pd,1);t3 = zeros(p.pd,1);t4 = zeros(odim,1); 
            t1(1:k) = y(1:k);
            for i =(k+1):p.pd
                t1(i) = s_linear(y(i),0.35);
            end
            t2(1:k) = t1(1:k);
            for i =(k+1):p.pd
                t2(i) = b_flat(t1(i),0.8,0.75,0.85);
            end
            for i =1:p.pd
                t3(i) = b_poly(t2(i),0.02);
            end
            for i =1:(odim-1)
                t4(i) = r_sum(t3((i-1)*k/(odim-1)+(1:(odim-1))'),2*(i-1)*k/(odim-1)+2*(1:(odim-1))');
            end
            t4(odim) = r_sum(t3((k+1):p.pd),2*((k+1):p.pd)');
            x = zeros(odim,1);
            for i = 1:(odim-1)
                x(i) = max(t4(odim),1)*(t4(i)-0.5) + 0.5;
            end
            x(odim) = t4(odim);                        
            S=2.*(1:odim)';
            D=1;
            hm=shape('convex',x,p.od);
            hm(odim)=mixedM(x,1,5);           
            f=D*x(odim)+S.*hm;
           % f=f';
        end
    end
    function hm=shape(name,x,M)
%         M=p.od;
        hm=zeros(M,1);
        switch name
            case 'linear'
                hm(1)=prod(x(1:(M-1)));
                for m=2:(M-1)
                    hm(m)=prod(x(1:(M-m)))*(1-x(M-m+1));
                end
                hm(M)=1-x(1);
            case 'convex'
                hm(1)=prod(1-cos(pi*x(1:(M-1))./2));
                for m=2:(M-1)
                    hm(m)=prod(1-cos(pi*x(1:(M-m))./2)).*(1-sin(pi*x(M-m+1)/2));
                end
                hm(M)=1-sin(pi*x(1)/2);
            case 'concave'
                 hm(1)=prod((sin(pi*x(1:(M-1))./2)));
                for m=2:(M-1)
                    hm(m)=prod((sin(pi*x(1:(M-m))./2)))*(cos(pi*x(M-m+1)/2));
                end
                hm(M)=cos(pi*x(1)/2);                               
        end                
    end
    function sp=mixedM(x,a,A)
        sp=(1-x(1)-(cos(2*A*pi*x(1)+pi/2))/(2*A*pi))^a;
    end
    function sp=discM(x,a,b,A)
        sp=1-(x(1)^a)*(cos(A*(x(1)^b)*pi))^2;
    end


    function t=b_poly(y,a)
        t=y.^a;
    end
    function t=b_flat(y,A,B,C)
        t=A+min(0,floor(y-B)).*((A.*(B-y))./B)-min(0,floor(C-y)).*(((1-A).*(y-C))./(1-C));
    end
    function t=b_param(y,y1,A,B,C)
        t = y^(B + (C-B) * v(y1));
        function vv=v(x)
            vv = A - (1-2*x)*abs(floor(0.5-x)+A);
        end        
    end
    function t=s_linear(y,A)
        t = abs(y-A)./(abs(floor(A-y)+A));
    end
    function t=s_decept(y,A,B,C)
        t = 1+(abs(y-A)-B).*((floor(y-A+B).*(1-C+(A-B)/B))./(A-B)+(floor(A+B-y).*(1-C+(1-A-B)/B))./(1-A-B)+1/B);
    end
    function t=s_multi(y,A,B,C)
        t = (1 + cos((4*A+2)*pi*(0.5-(abs(y-C))/(2*(floor(C-y)+C)))) +4*B*((abs(y-C))/(2*(floor(C-y)+C)))^2)/(B+2);
    end
    function t=r_sum(y,w)
        t=dot(y,w)/sum(w);
    end
    function t=r_nonsep(y,A)
        ny = length(y);
        tp1 = 0;
        for j =1:ny
            tp1 = tp1 + (y(j) + sum(abs(y(j) - y(mod(1+j + (0:(A-2)),ny)+1))));
        end
        t = (tp1)/((ny/A)*ceil(A/2)*(1+2*A-2*ceil(A/2)));
    end