%BY-AMARTYA SV
%Q.A planar 2D nozzleis considered. The flow is steady and frictionless and
%the density of the fluid is constant.Use the backward staggered grid with
%n pressure nodes and (n-1) velocity nodes.The stagnation pressure is given
%at the inlet and the static pressure is specified at the exit.Using the
%SIMPLE algorithm write down discretised momentum and pressure correction
%equations and solve for the unknown pressures and velocities. Check wether
%continuity is satisfied.
close;
clear all;
n=input('enter the number of grid points you want the domain split into ')
relax=input('enter the relaxation factor to be use to under relax pressure and velocity at the end of each iteration ')
%density and length of nozzle
rho=1;
L=2;
%initializing pressure and velocity. The pressure values at n nodes is
%stored as a (n)*1vector and the velocities are stored as (n-1)*1 vector
press=zeros(n,1);
vel=zeros(n-1,1);
%grid spacing
delx=L/(n-1);
%cross-sectional areas at pressure nodes
Ap=zeros(n,1);
for i=1:n
    Ap(i)=(1-(0.4*(i-1)*delx))/2;
end    
Ap;
%cross-sectional areas at velocity nodes
Av=zeros(n-1,1);
 for i=1:n-1
     Av(i)=(1-(0.4*delx*(1/2+(i-1))))/2;
 end    
 Av;
%initializing velocity values at nodes using an assumed mass flow rate of
%1kg/s
init_vel=zeros(n-1,1);
m=1;
for i=1:n-1
    init_vel(i)=m/(rho*Av(i));
end    
init_vel;
%initializing pressure values at nodes assuming a linear variation from A
%to E
init_press=zeros(n,1);
for i=1:n
    init_press(i)=(20-(10*(i-1)*delx))/2;
end
init_press;
%initializing a vector to track convergance-momentum residual
%mr=zeros(q,1);
%to run the loop q times
k=1;
momentum_residual=1;
while momentum_residual>0.001
    %discretised momentum equation:
    %mv contains the required coeffecients of all nodes. Each row of mv
    %contains the 4 co-effecients ap,aw,ae,su for that particular node.The
    %vector D contains the paramater d defined as A/ap for use in the
    %pressure correction equation.
    mv=zeros(n-1,4);
    D=zeros(n-1,1);
    %the nodes are split into 3 cases:interior nodes,the first node and the
    %last node.
    for i=1:n-1
        if i==1
           
            Fw=rho*((init_vel(i)*Av(i))/Ap(i))*Ap(i);
            Fe=rho*((init_vel(i)+init_vel(i+1))/2)*Ap(i+1);
            aw=0;
            ae=0;
            ap=Fe+(Fw*0.5*((Av(i)/Ap(i))^2));
            d=Av(i)/ap;
            D(i)=d;
            Su=(10-init_press(i+1))*Av(i)+(Fw*(Av(i)/Ap(i))*init_vel(i));
            mv(i,:)=[ap,aw,ae,Su];

        elseif i==n-1
            Fw=rho*((init_vel(i)+init_vel(i-1))/2)*Ap(i);
            Fe=m;
            aw=Fw;
            ae=0;
            ap=aw+ae+(Fe-Fw);
            Su=(init_press(i)-init_press(i+1))*Av(i);
            d=Av(i)/ap;
            D(i)=d;
            mv(i,:)=[ap,aw,ae,Su];

        else
            Fw=rho*((init_vel(i)+init_vel(i-1))/2)*Ap(i);
            Fe=rho*((init_vel(i)+init_vel(i+1))/2)*Ap(i+1);
            aw=Fw;
            ae=0;
            ap=aw+ae+(Fe-Fw);
            Su=(init_press(i)-init_press(i+1))*Av(i);
            d=Av(i)/ap;
            D(i)=d;
            mv(i,:)=[ap,aw,ae,Su];
        end
    end    
    mv;
    D;
    %The system of equation is solved as a linear system AX=b.X is the
    %velocities at the staggered nodes.mv contains the information for A
    %and b.mv is rearranged to yield a set of equations
    h=zeros(n-1);
    b=zeros(n-1,1);
    for i=1:n-1
        if i==1
            h(1,1)=mv(1,1);
            b(1)=mv(1,4);
        else
            h(i,i)=mv(i,1);
            h(i,i-1)=mv(i,2)*(-1);
            b(i)=mv(i,4);
        end
    end
    h;
    %calculating the momentum residual-The momentum residual is the
    %imbalance between the RHS and LHS of the discretized momentum
    %equ.Ideally as the iterations increase the residual has to tend to
    %zero.I've computed it here by substituting the initial velocity as
    %u.The final value for momentum residual is obtained by summing the
    %absolute values across all nodes.
    mom_res=h*init_vel-b;
    mom_res=abs(mom_res);
    momentum_residual=sum(mom_res);
    mr(k)=momentum_residual;
    
    b;
    vel_star=zeros(n-1,1);
    vel_star=pinv(h)*b;
    %thus the velocity is obtained from the discretised momentum
    %equation.These velocities correspond to v*
    %Now the pressure correction equation has to be solved.mp and B contain the
    %required coeffecients
    mp=zeros(n,3);
    B=zeros(n,1);
    for i=1:n
        if i==1
            mp(i,:)=[0,0,0];
            B(i)=0;
        elseif i==n
            mp(i,:)=[0,0,0];
            B(i)=0;
        else
            aw=rho*D(i-1)*Av(i-1);
            ae=rho*D(i)*Av(i);
            Fw=rho*vel_star(i-1)*Av(i-1);
            Fe=rho*vel_star(i)*Av(i);
            ap=aw+ae;
            b_=Fw-Fe;
            mp(i,:)=[ap,aw,ae];
            B(i)=b_;
        end
    end    
    mp;
    B;
    %solving the linear system as before to obtain the pressure corrections
    J=zeros(n);
    for i=1:n
        if i==1
            J(1,1)=1;
        elseif i==n
            J(n,n)=1;
        else
            J(i,i)=mp(i,1);
            J(i,i-1)=mp(i,2)*(-1);
            J(i,i+1)=mp(i,3)*(-1);

        end
    end
    press_corr=pinv(J)*B;
    %correcting the nodal pressures
    press=init_press+press_corr;
    %correcting velocities 
    for i=1:n-1
        vel(i)=vel_star(i)+D(i)*(press_corr(i)-press_corr(i+1));
    end
    vel;
    %corrected nodal pressure at A
    press(1)=10-(0.5*rho*(vel(1)^2)*((Av(1)/Ap(1))^2));
    press;
    
   
    %initialising velocity and pressure for next iteration with an
    %under_relaxation factor of 0.8 for both
    init_vel=(1-relax)*init_vel+(relax*vel);
    init_press=(1-relax)*init_press+(relax*press);
    %calculating the mass flow rate to cross check with the answer obtained
    %from bernoullis equation
    m=sum(rho*(init_vel.*Av))/(n-1);
    
    k=k+1;
end
%Table of results. Pressure displayed in first column,velocity in second
%and mass flow rate in 3rd.The mass flow rate is same at all nodes and
%hence displayed only once in the first row.The rest of the rows in the
%column are set to zero
Data=zeros(n,3);
Data(:,1)=press;
Data(n,1)=0;
Data(1:n-1,2)=vel;
Data(1,3)=m;
varnames={'pressure_values','velocity_values','mass_flow_rate'};
T = table(Data(:,1),Data(:,2),Data(:,3), 'VariableNames',varnames)
%To create the momentum residual plot
No_of_iterations_taken_to_converge=length(mr)
X=zeros(No_of_iterations_taken_to_converge,1);
for i=1:No_of_iterations_taken_to_converge
    X(i)=i;
end
mr;
plot(X,mr,'r-')
xlabel('iterations')
ylabel('momentum residual')



        

        
    

    
        
   
    
    
    
    
    
    
    
    
    
    
    
    

   
    
   
