capacitancia_membrana = 1;           
potencial_reposo = -65;              
max_conductancia_Na = 120;           
max_conductancia_K = 36;             
max_conductancia_fuga = 0.3;        
potencial_equilibrio_Na = 50;        
potencial_equilibrio_K = -77;       
potencial_equilibrio_fuga = -54.387; 


corriente_externa = [zeros(1, 500), 10*ones(1, 1500), zeros(1, 3500)];  


potencial_membrana = potencial_reposo; 
activacion_Na = 0.05; 
inactivacion_Na = 0.6; 
activacion_K = 0.32;


tiempo_final = 50; 
paso_temporal = 0.01;
tiempo = 0:paso_temporal:tiempo_final;


potencial_membrana_registro = zeros(1, length(tiempo));


for i = 1:length(tiempo)
   
    corriente_Na = max_conductancia_Na * activacion_Na^3 * inactivacion_Na * (potencial_membrana - potencial_equilibrio_Na);
    corriente_K = max_conductancia_K * activacion_K^4 * (potencial_membrana - potencial_equilibrio_K);
    corriente_fuga = max_conductancia_fuga * (potencial_membrana - potencial_equilibrio_fuga);

    
    dVdt = (corriente_externa(i) - corriente_Na - corriente_K - corriente_fuga) / capacitancia_membrana;
    dn_dt = 0.01 * (potencial_membrana + 55) / (1 - exp(-(potencial_membrana + 55) / 10)) * (1 - activacion_K) - 0.125 * exp(-(potencial_membrana + 65) / 80) * activacion_K;
    dm_dt = 0.1 * (potencial_membrana + 40) / (1 - exp(-(potencial_membrana + 40) / 10)) * (1 - activacion_Na) - 4 * exp(-(potencial_membrana + 65) / 18) * activacion_Na;
    dh_dt = 0.07 * exp(-(potencial_membrana + 65) / 20) * (1 - inactivacion_Na) - 1 / (1 + exp(-(potencial_membrana + 35) / 10)) * inactivacion_Na;

    
    potencial_membrana = potencial_membrana + paso_temporal * dVdt;
    activacion_K = activacion_K + paso_temporal * dn_dt;
    activacion_Na = activacion_Na + paso_temporal * dm_dt;
    inactivacion_Na = inactivacion_Na + paso_temporal * dh_dt;

    
    potencial_membrana_registro(i) = potencial_membrana;
end


figure;
plot(tiempo, potencial_membrana_registro);
title('Modelo Hodgkin-Huxley');
xlabel('Tiempo (ms)');
ylabel('Potencial de membrana (mV)');

