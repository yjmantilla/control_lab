% Programa que obtiene y grafica en una sola gr�fica la respuesta temporal del sistema Simulink 'ejemplo2�	
% cuando la entrada toma valores entre 1.5 y 7.5, con incrementos unitarios, para observar el efecto de un  
% bloque de saturaci�n en un sistema de control. El programa calcula s�lo 100 puntos de simulaci�n, entre 0 y 
% 9.9 segundos,  los cuales se completan mediante un procedimiento de interpolaci�n para obtener unas gr�ficas 
% m�s suaves. El programa tambi�n calcula el valor de estado estable de cada respuesta y su sobre nivel % porcentual.

clear simout ts xs ys     % Borre una a una todas las variables del espacio de trabajo de ML . 
% No use los comandos clear all  o clear  sin nombre de variables
t=[0:0.1:9.9];    % Se define el tiempo de simulaci�n, equivalente exactamente a 100 instantes de simulaci�n    
% en un vector columna. 
k = 1.5; % Amplitud inicial de la se�al de entrada
ks = {};
numIter = 7;

for i=1:numIter
    ut = [t;k*ones(size(t))]';                        	% Se generan las se�ales tipo escal�n de entrada al sistema
    % notice ut should be : An array with time in the first column and the remaining columns each corresponding to an input port. 
    
    [tt,xx,yy]=sim('ejemplo2',9.9,[],ut); 	% Se corre el programa Simulink y los resultados se guardan en las variables tt,xx,yy.    			          
    yint=interp1(tt,yy,t-0.5);    % La interpolaci�n crea puntos adicionales para unos resultados m�s precisos
    y0(:,i)=yint;	      % Se generan las  siete respuestas  del sistema
    ks{i} = num2str(k);
    k=k+1;
end


sy=size(y0);				% Se halla la dimensi�n de la variable
t_stable = sy(1);
y0ss=y0(t_stable,:);  			% Valor de respuesta en estado estable
SP=100*(max(y0)-y0ss)./y0ss; 	 	% Sobre nivel porcentual. No pase por alto el punto
for i=1:numIter
    ks{i} = ['k:' ks{i} ' - ' 'SP:' num2str(SP(i))];
end
plot(t,y0,'LineWidth',2); grid on;	%Se grafica el resultado parcial y se le agrega malla a los ejes
legend(ks,'Location','southeast')