clc;
%     o=0;
    sn=0;
for m=340

%     o=o+1;
%     mrk={'-s','-^','-+','-*'}
    disp(m)
    berr=[];
%     sn=0;
    for sn=0:5:60
        snr=10.^(sn/10);
    BER=nchoosek((2*m-1),m)*(0.5/snr).^m;
    berr=[berr,BER];
    end
%     disp(berr)
    sn=0:5:60;
    hold on
    set(gca, 'YScale', 'log')
    plot(sn,berr,'--','LineWidth',1.2,'color','blue');
    ylabel('BER')
    xlabel('SNR(dB)')
    axis([3 inf 10.^(-10000000) 1])
    legend('M=340','Location','northeast')
    box on
    
end
% MT=1:M;
% plot(MT,berr);