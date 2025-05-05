%analisi 


x1=out.simout(:,1);
x2=out.simout(:,2);
u=out.simout(:,3);
error=out.simout(:,4);


% Parametri
M = size(x1, 1);  
tolerance = 0.05;           % Tolleranza (5%)
% Valore finale 
final_value = x1(M);  

settling_time = 0;  
for k = 2:M
    if abs(x1(k) - final_value) > tolerance * abs(final_value)
        settling_time = k;
    end
end

sottoelongazione=0;
if min(x1)<x1(M)
    sottoelongazione=abs(min(x1)-x1(M))*100;
end

fprintf('Tempo di assestamento: %.2f min\n', 15*settling_time/M);
%fprintf('Errore max: %.2f mmoli/L\n', max(abs(error(:)))*1000);
fprintf('Temperatura max: %.2f K\n', max(x2));
fprintf('Temperatura min: %.2f K\n', min(x2));
fprintf('Input max: %.2f K\n', max(u));
fprintf('Input min: %.2f K\n', min(u));
fprintf('sottoelongazione percentuale: %.2f %% \n', sottoelongazione);
fprintf('Errore regime: %.2f mmoli/L\n', max(error(M))*1000);
