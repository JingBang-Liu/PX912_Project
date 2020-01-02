data = importdata('data1.txt')
figure(1)
plot(data(:,1),data(:,2),'r')
hold on
plot(data(:,1),data(:,3),'b')
% hold on
% % plot(data(:,1),data(:,4),'g')
% legend('initial','end-1')
figure(2)
plot(data(:,1),1+0.4*cos(data(:,1)))