% Programa que obtiene y grafica en una sola gráfica la respuesta temporal del sistema Simulink 'ejemplo2’	
% cuando la entrada toma valores entre 1.5 y 7.5, con incrementos unitarios, para observar el efecto de un  
% bloque de saturación en un sistema de control. El programa calcula sólo 100 puntos de simulación, entre 0 y 
% 9.9 segundos,  los cuales se completan mediante un procedimiento de interpolación para obtener unas gráficas 
% más suaves. El programa también calcula el valor de estado estable de cada respuesta y su sobre nivel % porcentual.

clear simout ts xs ys     % Borre una a una todas las variables del espacio de trabajo de ML . 
% No use los comandos clear all  o clear  sin nombre de variables
t=[0:0.1:9.9];    % Se define el tiempo de simulación, equivalente exactamente a 100 instantes de simulación    
% en un vector columna. 
k = 1.5; % Amplitud inicial de la señal de entrada
ks = {};
numIter = 7;

for i=1:numIter
    ut = [t;k*ones(size(t))]';                        	% Se generan las señales tipo escalón de entrada al sistema
    % notice ut should be : An array with time in the first column and the remaining columns each corresponding to an input port. 
    
    [tt,xx,yy]=sim('ejemplo2',9.9,[],ut); 	% Se corre el programa Simulink y los resultados se guardan en las variables tt,xx,yy.    			          
    yint=interp1(tt,yy,t-0.5);    % La interpolación crea puntos adicionales para unos resultados más precisos
    y0(:,i)=yint;	      % Se generan las  siete respuestas  del sistema
    ks{i} = num2str(k);
    k=k+1;
end


sy=size(y0);				% Se halla la dimensión de la variable
t_stable = sy(1);
y0ss=y0(t_stable,:);  			% Valor de respuesta en estado estable
SP=100*(max(y0)-y0ss)./y0ss; 	 	% Sobre nivel porcentual. No pase por alto el punto
for i=1:numIter
    ks{i} = ['k:' ks{i} ' - ' 'SP:' num2str(SP(i))];
end
plot(t,y0,'LineWidth',2); grid on;	%Se grafica el resultado parcial y se le agrega malla a los ejes
legend(ks,'Location','southeast')