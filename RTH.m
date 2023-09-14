function [Cost, Pos, Convergence_curve]=RTH(N,Tmax,low,high,dim,fobj)
Xbestcost = inf; Xbestpos = rand(N,dim);
for i=1:N
    Xpos(i,:) = low+(high-low).*rand(1,dim);
    Xcost(i)=fobj(Xpos(i,:));
    if Xcost(i) < Xbestcost
        Xbestpos = Xpos(i,:);
        Xbestcost = Xcost(i);
    end
end
A=15;
R0=0.5;
r = 1.5;
for t=1:Tmax
    %%               1- High Soaring
    Xmean=mean(Xpos);
    TF = 1+sin(2.5-t/Tmax); % TF
    for i=1:N
        Xnewpos=Xbestpos+(Xmean-Xpos(i,:)).*Levy(dim)*TF;
        Xnewpos = max(Xnewpos, low);
        Xnewpos = min(Xnewpos, high);
        Xnewcost=fobj(Xnewpos);
        if Xnewcost<Xcost(i)
           Xpos(i,:) = Xnewpos;
           Xcost(i)= Xnewcost;
             if Xcost(i) < Xbestcost
              Xbestpos= Xpos(i,:);
             Xbestcost=Xcost(i); 
             end
        end
    end
    %%                2- Low Soaring
    Xmean=mean(Xpos);
    for i=1:N-1
        aa=randperm(N);
        Xpos=Xpos(aa,:);
        Xcost=Xcost(aa);
        [x y]=polr(A,R0,N,t,Tmax,r);
        StepSize=Xpos(i,:)-Xmean;
        Xnewpos = Xbestpos +(y(i)+x(i))*StepSize;
        Xnewpos = max(Xnewpos, low);
        Xnewpos = min(Xnewpos, high);
        Xnewcost=fobj(Xnewpos);
        if Xnewcost<Xcost(i)
           Xpos(i,:) = Xnewpos;
           Xcost(i)= Xnewcost;
           if Xcost(i) < Xbestcost
                     Xbestpos= Xpos(i,:);
                    Xbestcost=Xcost(i);
           end
        end
    end
    %%                3- Stopping & Swooping
    Xmean=mean(Xpos);
    TF = 1+0.5*sin(2.5-t/Tmax); 
    for i=1:N
        b=randperm(N);
        Xpos=Xpos(b,:);
        Xcost=Xcost(b);
        [x y]=polr(A,R0,N,t,Tmax,r);
        alpha=(sin(2.5-t/Tmax).^2);
        G=2*(1-(t/Tmax));
        StepSize1 = 1.*Xpos(i,:) - TF*Xmean;
        StepSize2= G.*Xpos(i,:)-TF*Xbestpos;
        Xnewpos = alpha*Xbestpos+x(i)*StepSize1+y(i)*StepSize2;
        Xnewpos = max(Xnewpos, low);
        Xnewpos = min(Xnewpos, high);
        Xnewcost=fobj(Xnewpos);
        if Xnewcost<Xcost(i)
           Xpos(i,:) = Xnewpos;
           Xcost(i)= Xnewcost;
           if Xcost(i) < Xbestcost
                     Xbestpos= Xpos(i,:);
                    Xbestcost=Xcost(i); 
           end
        end
    end
    Convergence_curve(t)=Xbestcost;
    Cost = Xbestcost;
    Pos = Xbestpos;

end

 function [xR yR]=polr(A,R0,N,t,MaxIt,r)
%// Set parameters
th = (1+t/MaxIt)*A*pi*rand(N,1);
R  =(r-t/MaxIt)*R0*rand(N,1); 
xR = R.*sin(th);
yR = R.*cos(th);
 xR=xR/max(abs(xR));
 yR=yR/max(abs(yR));
 
 function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;