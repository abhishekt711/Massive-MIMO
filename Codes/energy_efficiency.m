clc;
M=input('enter the no. of antennas at the base station');
K=input('enter the no of antenas at the user');

    %4dB
a=0.6310;
b=1.5849;
pd = makedist('Normal','mu',0,'sigma',0.5);
pd_ = makedist('Normal','mu',1,'sigma',0.5);
tm = truncate(pd,-pi/9,pi/9);
t_m = truncate(pd,a,b);


P=0.005;
pds=0.005;
L=K;
vnoise=0.01;
dkkk=.36; dk=.51;
ak=[];
rt=0;
thetak=0.09;
S=sin(thetak);
lambda=1;
ts=0.5;
rate_=[];
rate_r=[];
rt=0;
rtr=0;
Psn=0;
Psnr=0;
psig=0.5;
eff=[];
eff_r=[];
for m=1:M
    it=1;
count = 0;


    %4dB
a=0.6310;
b=1.5849;
pd = makedist('Normal','mu',0,'sigma',0.5);
pd_ = makedist('Normal','mu',1,'sigma',0.5);
t = truncate(pd,-pi/9,pi/9);
t_ = truncate(pd,a,b);
r = random(t,m,1);
r_ = random(t_,m,1);



hbti=transpose(r_.*exp(1j.*r));
z=(sqrt(1/2).*(rand(1,m)+1i*rand(1,m)));
n=(sqrt(1/2).*(rand(m,K)+1i*rand(m,K)));
zk=transpose(z);
v=var(n(:));

    P=-1j*(m-1)*(pi)*S;
    akm=(exp(P));
    ak=[ak;akm];
    P=0.05;
    
    hh=((dkkk.*ak)+(dk.*zk));
hk=diag(hh);
H=diag(hbti*hk);
H_=diag(hh);

    phi=eye([K,K]);
phik=phi(:,K);
L=K;
vnoise=0.01;
hkhat=(dkkk.*ak)+((dk.^2)*sqrt(L*P)./((dk.^2).*L*P+vnoise)).*(((dk*sqrt(L*P)).*zk)+n*phik);
% disp(hkhat);
clear mean
expectationsquare=[sqrt(mean(hkhat.^2))];
%disp('values of wk are');
wk=[(hkhat/expectationsquare)];
Ps_=(norm(sqrt(pds).*lambda.*transpose(hk).*H*wk)).^2;
Ps_r=(norm(sqrt(pds).*lambda.*transpose(hk).*H_*wk)).^2;
Psn=Psn+Ps_;
Psnr=Psnr+Ps_r;
Ps_S=(norm(sqrt(psig).*lambda.*transpose(hk).*H*wk)).^2;
Ps_Sr=(norm(sqrt(psig).*lambda.*transpose(hk).*H_*wk)).^2;
% disp(Psn)
% disp(Ps_S)
SINRr=(Ps_Sr/(Psnr+v));
SINR=(Ps_S/(Psn+v));
rt_=(1-ts)*log(1+SINR)/log(2);
rt_r=(1-ts)*log(1+SINRr)/log(2);
rtr=rtr+rt_r;
rt=rt+rt_;
rate_=[rate_,rt];
rate_r=[rate_r,rtr];
eff=[eff,rt/(15+(0.05*m))];
end
MT=1:M;
% figure(1)
% plot(MT(1:20:M),rate_r(1:20:M),'-o',MT(1:20:M),rate_(1:20:M),'-*','LineWidth',1.2,'color','black')
% legend('Without reciprocity error ','with reciprocity error','Location','northwest')
% 
% xlabel('no. of antenna');
% ylabel('Achievable Rate (bits/s/Hz)');
% % title('Achievable Rate');
% 
% figure(2)
% plot(MT(1:10:M),rate_r(1:10:M),'-o',MT(1:10:M),rate_(1:10:M),'-*','LineWidth',1.2,'color','black')
% legend('Without reciprocity error ','with reciprocity error','Location','northwest')
% axis([150 300 7 12]);
% xlabel('no. of antenna');
% ylabel('Achievable Rate (bits/s/Hz)');
% title('Achievable Rate');

figure(1)
semilogx(MT(1:1:M),10*eff(1:1:M),'b--','LineWidth',1.2,'color','blue')
% legend('Without reciprocity error ','Location','northwest')
% axis([150 300 7 12]);
xlabel('Number of Antennas at BS (M)');
ylabel('Energy Efficiency(bits/sec/joule)');
%title('Energy Efficiency optimal curve of massive mimo');






