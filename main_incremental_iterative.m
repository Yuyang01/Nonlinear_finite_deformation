close all
clear all
global mod1 mesh1 load1 el1 undeformed1

% 0: upsetting of a block, dead load
% 1: upsetting of a block, imposed displacements
% 2: compression of a slender beam, imposed displacements
% 3: compression of a slender beam, dead load
% 4: arch, dead load at center of the arch
% 44: arch, imposed displacement at center of the arch
% 5: arch, dead load near the supports
example=1;
material=1;
[dof_force, dof_disp, lambda, x_eq, CC0, CC1, force, codeLoad]=preprocessing(example,material);

%Equilibrate
options.n_iter_max=80;
options.tol_x=1.e-6;
options.tol_f=1.e-6;
options.info=3;
options.method=0; %0: vanilla Newton-Rapshon, 1: Newton-Rapshon
options.linesearch=1; % 0: off, 1: on. For method 3, automatically on.

% Options for Line Search
options.n_iter_max_LS=30; % Maximum number of iterations for Line Search
options.type_LS=1; % 1: Backtracking, 2: Matlab
options.TolX=1.e-4;     % For Matlab line search
options.alfa=0.3; % For backtracking
options.beta = .8; % For backtracking

%Setup the undeformed configuration
precompute;

history_E=[];
history_x=[];
history_delta=[];
history_F=[];
%loop on the load increments
for iload=1:length(lambda)
    %Define the boundary conditions
    x=x_eq;
    load1.force = force*lambda(iload);
    switch example
        case {1, 2}
            x(1:2:end)=x_eq(1:2:end)*lambda(iload)/lambda(max(iload-1,1));
            load1.fixedvalues = x(load1.dofCC);
        case 44
            load1.fixedvalues(end-3:end)=load1.fixedvalues0(end-3:end)+lambda(iload)*load1.disp_max;
    end
    
    %x=x+rand(size(x))*.001; %random perturbations
    
    %Solve the equilibrium nonlinear system of equations
    [x_eq,iflag,iter,E_eq] = Equilibrate(x,options);
    [E_eq,grad_eq] = Energy(x_eq,2);
    history_E(iload)=E_eq;
    history_x(iload,:)=x_eq;
    switch example
        case {0, 1, 2, 3}
            history_delta(iload)=x_eq(2*CC1(1)-1)-mesh1.x0(2*CC1(1)-1);
            history_F(iload)=sum(grad_eq(2*CC0'-1)); %Reaction
        case {4, 5}
            history_delta(iload)=mean(x_eq(dof_disp)-mesh1.x0(dof_disp));
            history_F(iload)=mean(load1.force(dof_force)); %Reaction
        case {44}
            history_delta(iload)=mean(x_eq(dof_disp)-mesh1.x0(dof_disp));
            history_F(iload)=sum(grad_eq(load1.dofCC(end-3:end))); %Reaction
        otherwise
            disp('Case not implemented')
    end
    %Plot the deformed equilibrium configuration
    figure(1)
    clf
    DibujaMalla(mesh1.T,mesh1.x0,x_eq,'r',1)
end


%Plot the deformation vs. force
figure(3)
plot(-(history_delta),abs(history_F),'ro-')
xlabel('\delta')
ylabel('Force')

%Solve the same problem using the linear theory of elasticity
[u,Reaction,delta] = linear_elasticity(codeLoad,lambda(end),example,CC1,dof_force,dof_disp);
hold on
plot([0 -(delta)],[0,abs(Reaction)],'k-')
legend('Nonlinear elasticity','Linear elasticity')
figure(1)
hold on
DibujaMalla(mesh1.T,mesh1.x0,mesh1.x0+u,'k',1)
