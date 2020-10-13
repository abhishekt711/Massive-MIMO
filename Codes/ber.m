clc;
    o=0;
for m=2:2:8

    o=o+1;
    mrk={'s--','^--','o--','*--'}
    disp(m)
    berr=[];
    sn=0;
    for sn=5:5:60
        snr=10.^(sn/10);
    BER=nchoosek((2*m-1),m)*(0.5/snr).^m;
    berr=[berr,BER];
    end
%     disp(berr)
    sn=5:5:60;
    hold on
    set(gca, 'YScale', 'log')
    plot(sn,berr,mrk{o},'LineWidth',1.2);
    ylabel('BER')
    xlabel('SNR(dB)')
    axis([5 inf 10.^(-10) 1])
    legend('M=2','M=4','M=6','M=8','Location','northeast')
    box on
    
end
% MT=1:M;
% plot(MT,berr);