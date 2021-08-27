num = xlsread('PI Control.xls')

figure (1)
plot(num(:,1),num(:,2),'LineWidth',2)
hold on
plot(num(:,3),num(:,4),'g','LineWidth',2)
plot(num(:,5),num(:,6),'r','LineWidth',2)
hold off
legend('Mild','Severe','Set Point')
title('PI Control')
xlabel('Time (s)')
ylabel('Temperature (C)')