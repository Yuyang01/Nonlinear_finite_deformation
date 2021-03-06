function [x_equil,iflag,iter,Ener] = Equilibrate(x,options)
global mod1 mesh1 load1 el1 undeformed1


err_plot=[];
err_plot1=[];
[x_short] = short(x);

switch options.method
  case 0,   %vanilla Newton-Raphson
    iter=0;
    err_x=100;
    err_f=100;
    [Ener,grad_E,Hess_E] = Ener_short(x_short,3);
    while (iter<=options.n_iter_max) & ...
        ( (err_x>options.tol_x) | ...
        (err_f>options.tol_f))
      iter=iter+1;
      dx = -Hess_E\grad_E;
      x_short=x_short+dx;
      [Ener,grad_E,Hess_E] = Ener_short(x_short,3);
      err_x=norm(dx)/norm(x_short);
      err_f=norm(grad_E);
      err_plot=[err_plot err_x];
      err_plot1=[err_plot1 err_f];
      %fprintf('Iteration %i, errors %e %e \n', iter,err_x,err_f)
    end
    %Check positive definiteness
    if options.info==3
      [V,D] = eig(Hess_E);
      D=diag(D);
      if ((min(D))<=-1e-6*abs(max(D)))
           fprintf('Warning, the Hessian has a negative eigenvalue \n')
           D=sort(D);
           disp(D(1:6))
      end
    end
  otherwise,
    error('This option does not exist');
end

[x_equil] = long(x_short);
if (iter<=options.n_iter_max)
  if options.info>0
    fprintf('The equilibration was successful in %i iterations \n',iter)
  end
  iflag=1;
else
  if options.info>0
    fprintf('The equilibration was not reached in %i iterations \n',iter)
  end
  iflag=0;
end

%Output
if options.info>=2
  figure(5)
  hold off
  semilogy([1:iter],err_plot,'ro-',[1:iter],err_plot1,'bo-');
  %pause
end


