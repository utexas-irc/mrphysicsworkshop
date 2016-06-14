function [M, theta, phi] = pulsesim(B1, t, Foff, pulseshape, animate)
%PULSESIM Simulate off-resonant RF Pulses
%
%M = pulsesim(B1, t, Foff, pulseshape, animate)
%
%B1         = Magnitude of B1, in uT (~16 for human MRI scanners, ~25 for animal
%             MRI scanners, ~1200 for NMR spectrometers are reasonable starting
%             points.
%
%t          = Length of the RF pulse, in seconds.
%
%Foff       = Off-resonant frequency.  The difference between the Larmor
%             frequency of the nuclei at the field chosen and that of the RF
%             pulse.
%
%pulseshape = The name of the shaped pulse to use for the demonstration.
%             hard = rectangular pulse
%             sinc3 = 3-lobe sinc (one signle-sided zero crossing)
%             sinc5 = 5-lobe sinc (two signle-sided zero crossings)
%             sinc7 = 7-lobe sinc (three signle-sided zero crossings)
%             gauss = Gaussian truncated at typical 0.01% max amplitude
%
%animate    =  0 => simply show the final result 
%              1 => animate the demo
%             -1 => plot nothing, only retun magetization vectors
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up some physical constants and initialize vectors %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tInc = 1e-7;
gamma = 267.5;     %(rad/s/uT)
tVec = tInc:tInc:t;
offResInc = Foff*2*pi*tInc;
M = zeros(3,numel(tVec));
M(3,1) = 1;
animStep = numel(tVec)/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the RF pulse shape and incremental flips %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(pulseshape)
    case('hard')
        flipInc = gamma*B1*t/numel(tVec);  % in rad
        for ii = 2:numel(tVec)
            M(:,ii) = rotM(M(:,ii-1)', flipInc, offResInc)';
        end
        
    case('sinc3')
        ttemp = -2*pi:4*pi/(numel(tVec)-1):2*pi;
        sinc3 = B1*sin(ttemp)./ttemp;
        clear ttemp;
        flipInc = gamma*sinc3*tInc;
        for ii = 2:numel(tVec)
            M(:,ii) = rotM(M(:,ii-1)', flipInc(ii), offResInc)';
        end
        
    case('sinc5')
        ttemp = -3*pi:6*pi/(numel(tVec)-1):3*pi;
        sinc5 = B1*sin(ttemp)./ttemp;
        clear ttemp;
        flipInc = gamma*sinc5*tInc;
        for ii = 2:numel(tVec)
            M(:,ii) = rotM(M(:,ii-1)', flipInc(ii), offResInc)';
        end
        
    case('sinc7')
        ttemp = -4*pi:8*pi/(numel(tVec)-1):4*pi;
        sinc7 = B1*sin(ttemp)./ttemp;
        clear ttemp;
        flipInc = gamma*sinc7*tInc;
        for ii = 2:numel(tVec)
            M(:,ii) = rotM(M(:,ii-1)', flipInc(ii), offResInc)';
        end
        
    case('gauss')
        gauss = B1*(gausswin(numel(tVec), 3.71));
                flipInc = gamma*gauss*tInc;
        for ii = 2:numel(tVec)
            M(:,ii) = rotM(M(:,ii-1)', flipInc(ii), offResInc)';
        end
        
    otherwise
        disp('Pulse shape not recognized');
end

theta = sprintf( '%0.2f', (atan2(norm(cross([0 0 1],M(:,end))),dot([0 0 1],M(:,end))))*180/pi);
phiVec = [M(1,end) M(2,end) 0];
phi = sprintf( '%0.2f', (atan2(norm(cross([0 1 0],phiVec)),dot([0 1 0],phiVec)))*180/pi);

%%%%%%%%%%%%%%%%%%%%
% Plot the results %
%%%%%%%%%%%%%%%%%%%%
if((animate ~= 1) && (animate ~= 0))
    return;
end
magPlot = figure;
[sx, sy, sz] = sphere(100);
objHand = surf(sx, sy, sz);
set(objHand, 'LineStyle', 'none', 'FaceColor', [0.784 0.816 0.831]);
hold on

% These are just the axis markers - x = red, y = green, z = blue
plot3([1 -1], [0 0], [0 0], 'Marker', '.', 'MarkerSize', 15, 'Color', [1 0 0]);
plot3([0 0], [1 -1], [0 0], 'Marker', '.', 'MarkerSize', 15, 'Color', [0 1 0]);
plot3([0 0], [0 0], [1 -1], 'Marker', '.', 'MarkerSize', 15, 'Color', [0 0 1]);
set(gca, 'XLim', [-1 1], 'YLim', [-1,1], 'ZLim', [-1 1]);
set(gca, 'CameraPosition', [9.4318 7.1850 12.6261]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Animate for the demo if requested %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(animate == 1)
    %animate final result
    aPlot = plot3(M(1,1), M(2,1), M(3,1), 'Linewidth', 3);
    pPlot = plot3(M(1,1), M(2,1), M(3,1), 'Marker', '.', 'MarkerSize', 15, 'Color', [0.8 0 0]);
    movegui('west');
    
    timePlot = figure;
    pos = get(timePlot, 'Position');
    set(timePlot, 'Position', [pos(1) pos(2) 900 pos(4)]);
    movegui('east');
    subplot(1,3,1);
    xMagLine = plot(tVec(1), M(1,1), 'LineWidth', 2.5);
    xlabel('Time (s)', 'FontSize', 14);
    title('X-Magnetization', 'FontSize', 16, 'FontWeight', 'Bold');
    set(gca, 'XLim', [0 t],'YLim', [-1 1]);
    subplot(1,3,2);
    yMagLine = plot(tVec(1), M(2,1), 'LineWidth', 2.5);
    xlabel('Time (s)', 'FontSize', 14);
    title('Y-Magnetization', 'FontSize', 16, 'FontWeight', 'Bold');
    set(gca, 'XLim', [0 t],'YLim', [-1 1]);
    subplot(1,3,3);
    zMagLine = plot(tVec(1), M(3,1), 'LineWidth', 2.5);
    xlabel('Time (s)', 'FontSize', 14);
    title('Z-Magnetization', 'FontSize', 16, 'FontWeight', 'Bold');
    set(gca, 'XLim', [0 t],'YLim', [-1 1]);
    drawnow;   
    for ii = 2:animStep:numel(tVec)
        %figure(magPlot);
        set(aPlot, 'XData',M(1,1:ii),'YData', M(2,1:ii),'ZData', M(3,1:ii));
        set(pPlot, 'XData',M(1,ii),'YData', M(2,ii),'ZData', M(3,ii));
        %figure(timePlot);
        set(xMagLine, 'XData', tVec(1:ii), 'YData', M(1,1:ii));
        set(yMagLine, 'XData', tVec(1:ii), 'YData', M(2,1:ii));
        set(zMagLine, 'XData', tVec(1:ii), 'YData', M(3,1:ii));
        drawnow;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Or, just plot the final result %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    plot3(M(1,:), M(2,:), M(3,:), 'Linewidth', 3);
    movegui('west');
    figure;
    subplot(1,3,1);
    plot(tVec, M(1,:), 'LineWidth', 2.5);
    xlabel('Time (s)', 'FontSize', 14);
    title('X-Magnetization', 'FontSize', 16, 'FontWeight', 'Bold');
    set(gca, 'XLim', [0 t],'YLim', [-1 1]);
    subplot(1,3,2);
    plot(tVec, M(2,:), 'LineWidth', 2.5);
    xlabel('Time (s)', 'FontSize', 14);
    title('Y-Magnetization', 'FontSize', 16, 'FontWeight', 'Bold');
    set(gca, 'XLim', [0 t],'YLim', [-1 1]);
    subplot(1,3,3);
    plot(tVec, M(3,:), 'LineWidth', 2.5);
    xlabel('Time (s)', 'FontSize', 14);
    title('Z-Magnetization', 'FontSize', 16, 'FontWeight', 'Bold');
    set(gca, 'XLim', [0 t],'YLim', [-1 1]);
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 900 pos(4)]);
    movegui('east');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Otherwise, don't plot anything %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Report the final flip angle and phase %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(findobj(magPlot))
    figure(magPlot);
    title(['Final \theta = ' theta '\circ, Final \phi = ' phi '\circ'], 'FontSize', 16, 'FontWeight', 'Bold');
end



    function Mrot = rotM(M0, rB1, rOff)
        Mrot = M0*[1           0          0          ; ...
                   0           cos(rB1)  -sin(rB1)   ; ...
                   0           sin(rB1)   cos(rB1) ] ;
               
        Mrot = Mrot*[cos(rOff)  -sin(rOff)  0          ; ...
                     sin(rOff)   cos(rOff)  0          ; ...
                     0           0          1        ] ;

    end

end