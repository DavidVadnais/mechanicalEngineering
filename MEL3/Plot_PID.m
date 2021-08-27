num1 = tdfread('PID - Sheet1.tsv')

num = [num1.X_Value,num1.Var2, num1.Untitled,num1.X_Value1,...
    num1.Untitled1,num1.X_Value2,num1.SP];
figure (1)
plot(num(1:100,2),num(1:100,3),'LineWidth',2)
hold on
plot(num(1:100,4),num(1:100,5),'g','LineWidth',2)
plot(num(1:100,6),num(1:100,7),'r','LineWidth',2)
hold off
legend('Mild','Severe','Set Point')
title('PID Control')
xlabel('Time (s)')
ylabel('Temperature (C)')