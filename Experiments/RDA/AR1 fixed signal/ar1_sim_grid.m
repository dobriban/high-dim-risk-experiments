%% this code is a wrapper to perform simulations with an AR1 model
% the signal strength is set directly

cd('C:\Git\high-dim-risk-experiments\Experiments\RDA\AR1 extensive');
addpath '..\..\..\Code\'

%% run AR(1) simulations for a grid of parameters

%define grids
rho = [0.1; 0.5; 0.9];
gamma = [0.5; 1; 2];
alpha = [0.5; 1; 2];

%% test case
% gamma = 1;
% rho = 0.9;
% alpha = 1;

%% loop over parameters
a = {'-','--',':','-.'};
n_lambda = 100;
n = 5*1e2;
for i=1:length(rho)
    for j=1:length(gamma)
        for k=1:length(alpha)
            %perform experiment
            [lambda,risk,lambda_th,risk_th,oracle_error_rate] = run_ar1_sim(rho(i),gamma(j),alpha(k),n_lambda,n);
            
          
            %save results
            figure, hold on
            %normalize =  @(lambda)  (2*((lambda.*(lambda+4)).^(1/4))./(lambda.^(1/2)+(lambda+4).^(1/2))).^2;
            normalize =  @(lambda)  sqrt(lambda);
            h = plot(normalize(lambda),risk);
            set(h,'Linewidth',3,'LineStyle',a{1})
            h = plot(normalize(lambda_th),risk_th);
            set(h,'Linewidth',3,'LineStyle',a{2})
            h = plot(normalize(lambda_th),oracle_error_rate*ones(length(lambda_th),1));
            set(h,'Linewidth',3,'LineStyle',a{3})
            xlabel('$\sqrt{\lambda}$','Interpreter','LaTex')
            xlim([min(normalize(lambda)) max(normalize(lambda))]);
            ylim([0, 0.5]);
            ylabel('test error');
            %legend('Empirical P(err)','Theoretical P(err)','Oracle','Location','Best');
            set(gca,'fontsize',20)
            %titlestr = sprintf( 'rho=%.1f gamma=%.1f alpha=%.1f',rho(i),gamma(j),alpha(k));
            %title(titlestr);
            filename = sprintf( './ar1_rho_%.1f_gamma_%.1f_alpha_%.1f.pdf',rho(i),gamma(j),alpha(k));
            save2pdf(filename)
            %saveas(gcf, filename,'pdf');
            fprintf(['Saved Results to ' filename '\n']);
              
        end
    end
end