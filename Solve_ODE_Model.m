tf=360; 
Idodof=[1 5000 0 0 180 0 0 0 30 0 0 0 0 0 0 0 0]; 
vec = zeros(1,17);
for i = 1:17
    vec(i) = 1e-1;
end
options=odeset('RelTol',1e-1,'AbsTol',vec); 
[T,X]=ode45(@HCVIntraModel,[0 tf],Idodof,options); 
set(gcf,'position',[100 100 550 500])
figure(1)
set(0,'DefaultAxesFontSize',18)
RNA = X(:,1)+X(:,3)+X(:,11)+X(:,12)+2*X(:,13)+2*X(:,15)+X(:,16);
semilogy(T/24,RNA,'r','linewidth',3)
xlim([0 15])
xticks([0 5 10 15])
ylim([1 1e+8])
yticks([1 1e+2 1e+4 1e+6 1e+8])
hold on
semilogy(T/24,X(:,5),'b','linewidth',2.5)
hold on
semilogy(T/24,X(:,6),'color','0 0.5 0','linewidth',2.5)
hold on
semilogy(T/24,X(:,17),'k','linewidth',2.5)
set(figure(1),'Renderer','painters')
