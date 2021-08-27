num = xlsread('P Only Control.xls')

figure (1)
plot(num(1:100,2),num(1:100,3),'LineWidth',2)
hold on
plot(num(1:100,4),num(1:100,5),'g','LineWidth',2)
plot(num(1:100,6),num(1:100,7),'r','LineWidth',2)
hold off
legend('Mild','Severe','Set Point')
title('P Only Control')
xlabel('Time (s)')
ylabel('Temperature (C)')