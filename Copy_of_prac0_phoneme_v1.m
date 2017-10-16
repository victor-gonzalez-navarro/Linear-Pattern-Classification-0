% PRAC0
%  MC Feb. 2017
% Lectura de la base de datos phoneme
% Análisis de Gaussianidad: Gráfico y Análitico

%% OPCIONES
clear
close all
i_dib=0;                        %0 NO /1 SI: Dibuja señales
i_histfit=1;                    %0 NO /1 SI: HISTOGRAMAS
i_cdfplot=0;				    %0 NO /1 SI: calcula cdf
i_scplot=0;					    %0 NO /1 SI: SCATTERPLOT DE CARACTERISTICAS
i_kurskew=0;				    %0 NO /1 SI: calcula Kurtosis Y Skewness
i_plotnorm=0;				    %0 NO /1 SI: calcula PLOTNORM
%% Parámetros de la BD                 
% class name: Labels:
clases=['aa';'ao';'dc';'iy';'sh'];  % 1(aa);2(ao);3(dcl);4(iy);5(sh);
N_clas=5;
N_fft=256;						%256 (8KHz) 128 (4KHz), 64 (2KHz), 32(1khZ)
%% Database load
load BD_phoneme

%% Spectrum plot
if i_dib==1
    Frec_max=8*N_fft/256;			%Max frequency in KHz
    eje_frec=(0:N_fft-1)*Frec_max/N_fft;
    figure('name','LOG(Espectrum)')
    for i_clas=1:N_clas
        subplot(3,2,i_clas)
        hold on
        index=find(Labels==i_clas);
        for i1=1:length(index)
            plot(eje_frec,X(index(i1),1:N_fft));
        end
        hold off
        grid
        zoom on
        xlabel('frec(KHz)')
        ylabel(clases(i_clas,:));
    end
    subplot(3,2,N_clas+1)
    hold on
    i_color=['b' 'r' 'g' 'k' 'y'];
    for i_clas=1:N_clas
        index= Labels==i_clas;
        aux=mean(X(index,1:N_fft));
        plot(aux,i_color(i_clas));
    end
    hold off
    grid
    zoom on
    xlabel('Feature Number')
    ylabel('log espectro')
    title('Average');
    clear index aux i_color i_clas eje_frec Frec_max
end
clear i_dib N_fft

%% Feature selection: 6 características para analizar
V_coor=[25 50 75 100 125 150];
N_feat=length(V_coor);
if V_coor(1)~=0
    X=X(:,V_coor);  % Feature selection
end
%%  KURTOSIS-SKEWNESS
if i_kurskew==1
    figure('name','Kurtosis,Skewness')
    for i_clas=1:N_clas
        index=find(Labels==i_clas);
        subplot(2,5,i_clas)
        bar(V_coor,kurtosis(X(index,:))-3)
        ylabel('KURTOSIS')
        xlabel('feature number')
        title(clases(i_clas,:));
        grid
        subplot(2,5,N_clas+i_clas)
        bar(V_coor,skewness(X(index,:)))
        ylabel('SKEWNESS')
        xlabel('feature number')
        title(clases(i_clas,:));
        grid
    end
end

%% SCATTER PLOT
if i_scplot==1
    varNames = {'feat 1' 'feat 2' 'feat 3' 'feat 4' 'feat 5' 'feat 6'};
    figure('name','Scatter Plot')
    gplotmatrix(X,X,Labels,'bgrky',[],[],'on','hist',varNames,varNames)
    zoom on
end

%% HISTOGRAMA
i_class=3 % Clase que se analiza: 1(aa);2(ao);3(dcl);4(iy);5(sh);
if i_histfit==1
    figure('name','Histograma de la clase seleccionda')
    index=find(Labels==i_class);  %Seleccionar clase a analizar
     for i_feat=1:N_feat
        subplot(3,2,i_feat)
        histfit(X(index,i_feat))
        grid
        zoom on
        title(i_feat)
    end
end


length(index) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mediap = zeros(1,i_feat);
for j=1:i_feat
    for i=1:(length(index))
        mediap(j) = mediap(j) + X(index(i),j);
    end
end
media = mediap/(length(index));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cuant = zeros(1,i_feat);
for j=1:i_feat
    for i=1:(length(index))
        cuant(j) = cuant(j) + (X(index(i),j)- media(j))*(X(index(i),j)- media(j));
    end
end
s2 = cuant/(length(index)-1);
s = sqrt(s2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.05;
p = 1 - alpha/2;
df = length(index) - 1;
t = tinv(p,df);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intervalo1 = zeros(1,i_feat);
for j=1:i_feat
    intervalo1(j) = (media(j) - t*(s(j)/sqrt(length(index))));
end
intervalo1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intervalo2 = zeros(1,i_feat);
for j=1:i_feat
    intervalo2(j) = (media(j) + t*(s(j)/sqrt(length(index))));
end
intervalo2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_class
parameter = 1

x = zeros(1,length(index));
for i=1:(length(index))
        x(i) = X(index(i),parameter);
end

%figure(3), histfit(X(index,1));
%figure(2), histogram(x,30);
%[h,p,stats] = chi2gof(x,'Alpha',0.001,'NBins', 30, 'Emin', 0)
[h,p,stats] = chi2gof(x,'Alpha',0.001)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cdf
if i_cdfplot==1
    figure('name','cdfplot de la clase seleccionda')
    index=find(Labels==i_class);
    for i_feat=1:N_feat
        subplot(3,2,i_feat)
        Aux=X(index,i_feat);
        Yaux=linspace(min(Aux),max(Aux),500);
        plot(Yaux,cdf('Normal',Yaux,mean(Aux),std(Aux)),'r');
        hold on
        cdfplot(X(index,i_feat))
        title(i_feat)
        zoom on
    end
    clear Aux Yaux
end

%%  PLOTNORM
if i_plotnorm==1
    figure('name','Plotnorm de la clase seleccionda')
    index=find(Labels==i_class);
    for i_feat=1:N_feat
        subplot(3,2,i_feat)
        qqplot(X(index,i_feat))
        grid
        title(i_feat)
        zoom on
    end
end


