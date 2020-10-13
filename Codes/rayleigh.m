pd = makedist('Rician','s',0,'sigma',2)
pd1 = makedist('Rayleigh','b',2.2)
pd2 = makedist('Normal','mu',2,'sigma',1.5);
x = -10:.1:10;

plot(x,pdf(pd,x),x,pdf(pd1,x),x,pdf(pd2,x),'LineWidth',1.2);
legend('Rician fading','Rayleigh fading','Gaussian distribution','Location','northwest')