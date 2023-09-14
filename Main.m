%%  Red-Tailed Hawk Algorithm (RTH) code:
%   Main paper: Red-tailed hawk algorithm for numerical optimization and real-world problems
%   DOI: 10.1038/s41598-023-38778-3
%% Main code
clear all ; clc
Npop=30;                
Function_name='F10';     % Name of the test function that can be from F1 to F23 ( 
Tmax=1000;              
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
[Best_fit,Best_pos,Convergence_curve]=RTH(Npop,Tmax,lb,ub,dim,fobj);
figure('Position',[500 500 660 290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])
%Draw objective space
subplot(1,2,2);
semilogy(Convergence_curve,'Color','r')
title('Search space')
xlabel('Iteration');
ylabel('Best score obtained so far');
axis tight
grid on
box on
legend('RTH')
display(['The best solution is ', num2str(Best_pos)]);
display(['The best fitness value is ', num2str(Best_fit)]);