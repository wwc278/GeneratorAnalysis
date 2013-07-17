function [PtoHAggregate PtoHCorrelate DataOut n M] = ...
    gen_data_analysis(design_torque, max_heat, T)
%Generator Data Analysis
%   This file analyzes the maxon motor specifications found in variable
%   MotorData and outputs useful information in DataOut. Efficiency over a
%   range of torque and speeds, friction losses, heat losses, electrical
%   output power, etc.

clc
close all
load('MaxonMotorDataOG.mat','MotorData')
[row,col]=size(MotorData);


if ~exist('rpm_limit') || isempty(rpm_limit)
    rpm_limit=80000;
end
if ~exist('design_torque') || isempty(design_torque)
    design_torque=23;
end
if ~exist('max_heat') || isempty(max_heat)
    max_heat=6;         %Maximum heat output in Watts
end
if ~exist('T') || isempty(T)
    T=155;         %Maximum heat output in Watts
end



for motornumber=1:col
%for motornumber=[26]
    
    a=MotorData(:,motornumber);
    %calculate efficiency and torque
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % automatic input of constants----------------------------
    k_M=a{13};          %(mNm/A) torque constant
    k_n=a{14};          %(rpm/mNm) speed constant
    I_o=a{4}*1e-3;      %(A) no load current
    n_o=a{3};           %(rpm) no load speed
    gr=a{15};           %(rpm/mNm) speed/torque gradient
    U=a{2};             %(V) voltage input
    M_H=a{8};           %(mNm) stall torque
    R25=a{11};          %(ohms) terminal resistance phase to phase
    alpha_Cu=0.0039;    %(1/K) temperature coefficient of copper at 20degC
    M_N=a{6};           %(mNm) nominal torque
    %---------------------------------------------------------
    
    %fundamental equations for efficiency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %M_R=k_M*I_o;        %(mNm) friction torque (torque at no load current)
    
    %[n,M]=meshgrid(0:15000/100:15000,0:M_N/200:M_N);
    
    [n,M]=meshgrid(0:rpm_limit/100:rpm_limit,...
        0:100/200:100);
    I=M/k_M;            %current related to torque
    
    %find friction torque with new formula for speeds greater than the no
    %load speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x,y]=size(n);
    
    for i=1:x
        for j=1:y
            M_R1=k_M*I_o;
            M_R2=k_M^2*pi/R25/30000*(n(i,j)-n_o);
   %         if M_R1>M_R2
                M_R(i,j)=M_R1;
   %         else
   %             M_R(i,j)=M_R2;
    %        end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    R_T=R25*(1+alpha_Cu*(T-25));
    
    P_mech=(pi/3e4)*n.*M;
    P_fric=(pi/3e4)*n.*M_R;
    P_Jo=R_T.*I.^2;
    
    P_el=n/k_n.*I;
    
    eta=(P_mech-P_fric-P_Jo)./P_mech;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %compute power to heat ratio for current motor and compare the values
    %to other motors over a range of speeds and torques%%%%%%%%%%%%%%%%%%%%
    
    Heat=P_fric+P_Jo;
    PowerToHeat=(P_mech-P_fric-P_Jo)./(P_fric+P_Jo);
    
    if ~exist('PtoHAggregate')
        PtoHAggregate=zeros(x,y);
    end
    if ~exist('PtoHCorrelate')
        PtoHCorrelate=cell(x,y);
    end
    
    for i=1:x
        for j=1:y
            if Heat(i,j)<max_heat && PowerToHeat(i,j)>PtoHAggregate(i,j)
                PtoHAggregate(i,j)=PowerToHeat(i,j);
                PtoHCorrelate{i,j}=a(1);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %change all negative efficiencies to NaN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [row,col]=size(eta);
    for j=1:col
        for k=1:row
            if eta(k,j)<0
                eta(k,j)=NaN;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %find max efficiency at each torque and corresponding values%%%%%%%%%%%
    for j=1:row
        [pks,ind]=max(eta(j,:));
        if isempty(pks)
            eta_output(j)=nan;
            M_output(j)=nan;
            n_output(j)=nan;
            P_M(j)=nan;
            P_F(j)=nan;
            P_J(j)=nan;
            P_el(j)=nan;
        else
            eta_output(j)=pks;
            M_output(j)=M(j,ind);
            n_output(j)=n(j,ind);
            P_M(j)=P_mech(j,ind);
            P_F(j)=P_fric(j,ind);
            P_J(j)=P_Jo(j,ind);
            P_E(j)=P_el(j,ind);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %plot surface and contour for visualization purposes%%%%%%%%%%%%%%%%%%%
    %
    %     if motornumber==11
    %         figure
    %
    %         subplot(2,2,1)
    %         surf(M,n,eta,'EdgeColor','none')
    %         xlabel('Torque (mNm)')
    %         ylabel('Speed (RPM)')
    %         zlabel('Efficiency')
    %         caxis([0 1])
    %         grid on
    %         title([MotorData{1,motornumber} ' ' num2str(MotorData{2,motornumber}) 'V '...
    %             num2str(MotorData{18,motornumber}) 'W ' '@' num2str(T) 'degC'])
    %
    %
    %         subplot(2,2,3)
    %         [C,h]=contourf(M,n,eta);
    %         clabel(C,h)
    %         caxis([0 1])
    %         xlimits=get(gca,'Xlim');
    %         ylimits=get(gca,'Ylim');
    %         xtick = get(gca,'XTick');
    %         ytick = get(gca,'YTick');
    %         %title('Efficiencies at various Torques and Speeds')
    %         xlabel('Torque [mNm]')
    %         ylabel('Speed [rpm]')
    %         grid on
    %         ax1=gca;
    %         ax2 = axes('Position',get(gca,'Position'),...
    %             'XAxisLocation','top',...
    %             'YAxisLocation','right',...
    %             'Color','none',...
    %             'XColor','k','YColor','k');
    %         set(ax2,'Xlim',xlimits/a{13},'Ylim',ylimits/a{14})
    %         set(ax2,'XTick',[xtick/a{13}],...
    %             'YTick',[ytick/a{14}])
    %         xlabel('Current [A]')
    %         ylabel('Voltage [V]')
    %         axes(ax1)
    %
    %         subplot(2,2,2)
    %         plot(M_output,P_M-P_F-P_J)
    %         grid on
    %         xlabel('Torque [mNm]')
    %         ylabel('Power [W]')
    %         title('Power Out vs. Torque')
    %         if max(P_M-P_F-P_J)>0
    %             axis([0 M_output(end) 0 max(P_M-P_F-P_J)])
    %         end
    %         %     subplot(2,2,4)
    %         %     x=M_output(1:round(length(M_output)/10):end);
    %         %     y=[P_F(1:round(length(P_F)/10):end);...
    %         %         P_J(1:round(length(P_J)/10):end);...
    %         %         P_M(1:round(length(P_M)/10):end)-...
    %         %         P_F(1:round(length(P_F)/10):end)-...
    %         %         P_J(1:round(length(P_J)/10):end)]';
    %         %
    %         % %     P_M(1:round(length(P_M)/10):end)-...
    %         % %         P_F(1:round(length(P_F)/10):end)-...
    %         % %         P_J(1:round(length(P_J)/10):end)
    %         %
    %         %     bar(x(1:end-1),y(1:end-1,:))
    %         %     xlim([0 M_output(end)])
    %         %     grid on
    %         %     xlabel('Torque [mNm]')
    %         %     ylabel('Power Loss [W]')
    %         %     title('Power Losses vs. Torque')
    %         %     legend('Friction Losses','Heat Losses','Available Power')
    %
    %
    %         subplot(2,2,4)
    %         plot(M_output,P_F,M_output,P_J)
    %         xlim([0 M_output(end)])
    %         grid on
    %         xlabel('Torque [mNm]')
    %         ylabel('Power Loss [W]')
    %         title('Power Losses vs. Torque')
    %         legend('Friction Losses','Heat Losses')
    %
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %         figure
    % %         subplot(2,2,1)
    % %         plot(M_output,eta_output)
    % %         grid on
    % %         xlabel('Torque [mNm]')
    % %         ylabel('Efficiency')
    % %         title('Max. Efficiency possible vs. Torque')
    % %         axis([0 10 0 1])
    % %
    % %         subplot(2,2,2)
    % %         plot(M_output,P_F,M_output,P_J)
    % %         grid on
    % %         xlabel('Torque [mNm]')
    % %         ylabel('Power Loss [W]')
    % %         title('Power Losses vs. Torque')
    % %         xlim([0 10])
    % %         legend('Friction Losses','Heat Losses')
    % %
    % %         subplot(2,2,3)
    % %         plot(M_output,P_M-P_F-P_J)
    % %         grid on
    % %         xlabel('Torque [mNm]')
    % %         ylabel('Power [W]')
    % %         title('Power Out vs. Torque')
    % %         if max(P_M-P_F-P_J)>0
    % %             axis([0 10 0 max(P_M-P_F-P_J)])
    % %         end
    % %
    % %         ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1]...
    % %             ,'Box','off','Visible','off','Units','normalized', ...
    % %             'clipping' , 'off');
    % %
    % %         text(0.5, 1,[['\bf' MotorData{1,i} ' ' num2str(MotorData{2,i})]...
    % %             ['V ' num2str(MotorData{18,i}) 'W @' num2str(T) 'degC']],...
    % %             'HorizontalAlignment','center','VerticalAlignment', 'top')
    % %
    % %
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %store information into DataOut%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DataOut(1,1)={'Motor Name'};
    DataOut(1,2)={'Rated Voltage'};
    DataOut(1,3)={'Rated Watts'};
    DataOut(1,4)={'Heat Output'};
    DataOut(1,5)={'Power Output'};
    DataOut(1,6)={'Power to Heat Ratio'};
    DataOut(1,7)={'RPM'};
    %     DataOut(1,8)={''};
    %     DataOut(1,9)={''};
    %     DataOut(1,10)={''};
    %     DataOut(1,11)={''};
    %     DataOut(1,12)={''};
    %     DataOut(1,13)={''};
    %     DataOut(1,14)={''};
    
    
    [C,I]=min(abs(M_output-design_torque));
    J=1;
    
    for i=1:col
        if Heat(I,i)<max_heat && PowerToHeat(I,i)>0 &&...
                PowerToHeat(I,i)>=PowerToHeat(I,J)
            J=i;
        end
    end
    
    
    if J~=1
        [rows cols]=size(DataOut);
        DataOut(rows+1,1)=a(1);
        DataOut(rows+1,2)=a(2);
        DataOut(rows+1,3)=a(18);
        DataOut(rows+1,4)={Heat(I,J)};
        DataOut(rows+1,5)={Heat(I,J)*PowerToHeat(I,J)};
        DataOut(rows+1,6)={PowerToHeat(I,J)};
        DataOut(rows+1,7)={n(I,J)};
        %         DataOut(8,row+1)={P_J};
        %         DataOut(9,row+1)={P_M-P_F-P_J};
        %         DataOut(10,row+1)=a(13);
        %         DataOut(11,row+1)=a{14};
        %         DataOut(12,row+1)={M};
        %         DataOut(13,row+1)={n};
        %         DataOut(14,row+1)={eta};
    end
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



figure
contourf(M,n,PtoHAggregate,'EdgeColor','none')
colorbar
grid on
title({'Power to Heat Ratio over a span of Torques and Speeds';...
    ['Max. Heat = ' num2str(max_heat) ', Temp. = ' num2str(T) 'degC']})
xlabel('Torque [mNm]')
ylabel('Speed [rpm]')

DataOut    %displays DataOut for debugging

[C,I]=max(cell2mat(DataOut(2:end,6)));
if ~isempty(I)
    I=I+1;
    Optimal_generator=[DataOut{I,1} ' ' num2str(DataOut{I,2}) 'V ' ...
        num2str(DataOut{I,3}) 'W with a heat output of ' ...
        num2str(DataOut{I,4}) 'W and a power output of ' ...
        num2str(DataOut{I,5}) 'W and a power to heat ratio of '...
        num2str(DataOut{I,6})]
    Design_limits=[ num2str(design_torque) ...
        'mNm, ' num2str(max_heat) 'W of heat, ' num2str(T) 'degC']
else
    disp('No motor found! Try increasing the heat constraint.')
end



end

