function [] = plotsignal(time, signal)

    figure
    subplot(2,2,1)
    plot(time,abs(signal),"LineWidth",1.1);
    grid on
    grid minor
    hold on
    xlabel("time [s]");
    ylabel("signal [amplitude]");
    title("Signal in Time")
    subplot(2,2,3)
    plot(time,angle(signal),"LineWidth",1.1);
    grid on
    grid minor
    hold on
    xlabel("time [s]");
    ylabel("signal [phase]");

    %% Fourier Transform
    dt = time(2)-time(1);
    Nf = 2*length(time);

    [Signal_FT,freq] = dft(signal,time,Nf);
    Signal_FT = Signal_FT*dt;
    df = freq(2)-freq(1);

    subplot(2,2,2)
    plot(freq,abs(Signal_FT),"LineWidth",1.1);
    grid on
    grid minor
    hold on
    xlabel("frequency [Hz]");
    ylabel("Signal Transform [amplitude]");
    title("Fourier Transform")

    subplot(2,2,4)
    plot(freq,(angle(Signal_FT)),"LineWidth",1.1);
    grid on
    grid minor
    hold on
    xlabel("frequency [Hz]");
    ylabel("Signal Transform [phase]");


    %% test
    %figure
    %subplot(1,2,1)
    %plot(time,real(abs(signal).*exp(1i*angle(signal))),"LineWidth",1.1);
    %grid on
    %grid minor
    %hold on
    %    y = idft(Signal_FT,freq',time');
    %    y = y*df*Nf;
    %plot(time, real(y),"LineWidth",1.1);
    %xlabel("time [s]");
    %ylabel("signal [amplitude]");
    %title("Signal in Time");
    %legend("Initial Signal", "DFT->IDFT Signal")
%
    %subplot(1,2,2)
    %plot(freq,real(abs(Signal_FT).*exp(1i*angle(Signal_FT))),"LineWidth",1.1);
    %grid on
    %grid minor
    %hold on
    %xlabel("frequency [Hz]");
    %ylabel("Signal Transform [amplitude]");
    %title("Fourier Transform")