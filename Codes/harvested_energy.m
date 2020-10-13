clc;
% M=input('enter the no. of antennas at the base station');
% K=input('enter the no of antenas at the user');
M=400;
K=6;
mrk={'-s','-^','-o','-*'};
for o=1:4
    if o==1
% 1db
a=0.8913;
b=1.1220;
    end
     if o==2
%4dB
a=0.6310;
b=1.5849;
     end
    if o==3
%6db
a=0.5012;
b=1.9953;
    end
     if o==4
%8db
a=0.3981;
b=2.5119;
    end     
pd = makedist('Normal','mu',0,'sigma',0.5)
t = truncate(pd,-pi/9,pi/9);
% legend({'Normal','Truncated'},'Location','NE')
r = random(t,M,1);

pd_ = makedist('Normal','mu',1,'sigma',0.5)
t_ = truncate(pd,a,b);
% legend({'Normal','Truncated'},'Location','NE')
r_ = random(t_,M,1);

hbti=transpose(r_.*exp(1j.*r));










%%noise generation


it=10000;
count = 0;
nn=0;
zkt=0;
zzgg=0;
zgg=0;
for i=1:it
y1=randn(M,K);
y2=1i*randn(M,K);

nk=(sqrt(1/2)*(y1+y2));
nit=[nk];
zg=(sqrt(1/2).*(rand(1,M)+1i*rand(1,M)));
zzg=[zg];
zkf=transpose(zzg);
zkt=zkt+zkf;
nn=nn+nit;
zgg=zgg+zg;
zzgg=zzgg+zzg;
end
n=nn/it;
zk=zkt/it;
z=zgg/it;
zz=zzgg/it;

thetak=.09;
S=sin(thetak);
ak=[];
for m=1:M-1
    P=-1j*(m-1)*(pi)*S;
    ak(m)=(exp(P));
    ak=[ak;ak(m)];
end




dkkk=.36; dk=.51;
hh=((dkkk.*ak)+(dk.*zk));
hk=diag(hh);
H=diag(hbti).*hk;
% disp(H);

%generating PILOT SEQUENCES
phi=eye([K,K]);
phik=phi(:,K);

P=0.05;
L=K;
vnoise=0.01;
hkhat=(dkkk.*ak)+((dk.^2)*sqrt(L*P)./((dk.^2).*L*P+vnoise)).*(((dk*sqrt(L*P)).*zk)+n*phik);
%disp(hkhat);
clear mean
expectationsquare=[sqrt(mean(hkhat.^2))];
%disp('values of wk are');
wk=[(hkhat/expectationsquare)];
%disp(wk);
% sk=1;
% pk=1;
% yk=hH*wk;
hH_=ctranspose(hk);
hH=ctranspose(H);
Ek=0;
m=1;
E=[];
E_=[];
% n = input(' Enter n: ');
n=100000;
count = 0;
Qs=0;
Qs_=0;


ps=0.6;
E=[];
E_=[];
Qs=0;
Qs_=0;
m=1;
for i=1:n,
delk=0.4*(1-ps)*0.005;
    while m<M+1
    se=(wk);
    Ek=se+Ek;
    Qk=delk.*((hH*Ek).^2);
    pe=(sum(abs(Qk)));
    E=[E;pe];
    Qk_=delk.*((hH_*Ek).^2);
    pe_=(sum(abs(Qk_)));
    E_=[E_;pe_];
    m=m+1;
    if m > M
%         ps=ps-0.3;
%         m=1;
%         if ps<0
%             break;
%         end
        break;
    end
    end

    Qs=Qs+E;
    Qs_=Qs_+E_;


end

Qsk=Qs./n;

Qsk_=Qs_./n;


    

 MT=[1:1:M];


ptt=[0;diff(Qsk)];
ptt_=[0;diff(Qsk_)];
figure(1)

hold on;
axis([100 384 0 inf])
cc=ptt/M;
box on


plot( MT(1:10:M), cc(1:10:M),mrk{o},'LineWidth',1.2)

ylabel('Harvested Energy (J)')
xlabel('Number of Antennas at BS (M)');

% title('Comparison of channel reciprocity error')
% disp('ptt');
% disp(ptt/M);
%disp(diff(Qsk));
end
c=ptt_/M;
plot( MT(1:10:M), c(1:10:M),'-','LineWidth',1.2)
legend('1dB error','4dB error','6dB error','8dB error','without reciprocity error','Location','northwest')

