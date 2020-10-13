pd = makedist('Rician','s',0,'sigma',2)
pd1 = makedist('Rician','s',0,'sigma',3)
x = 0:.1:25;
plot(x,pdf(pd,x),x,pdf(pd1,x),'LineWidth',1.2);
legend('K=-\infty dB','K=6dB','Location','northeast')