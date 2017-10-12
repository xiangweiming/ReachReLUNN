clc
clear
%% Neural network
%Randomly generate neural network 
%numNueron = [3,5,5,5,5,2];
%network.activefun = {'poslin','poslin','poslin','poslin','purelin'};

%numNueron = [5,9,9,9,9,3];
%network.activefun = {'poslin','poslin','poslin','poslin','purelin'};
%network.activefun = {'purelin','purelin','purelin','purelin','purelin'};

numNueron = [3,7,7,7,7,7,7,7,2];
network.activefun = {'poslin','poslin','poslin','poslin','poslin','poslin','poslin','purelin'};

%numNueron = [3,9,9,9,9,9,9,9,9,9,2];
%network.activefun = {'poslin','poslin','poslin','poslin','poslin','poslin','poslin','poslin','poslin','purelin'};

numLayer = length(numNueron)-1;
for n = 1:1:numLayer
    W{n} = randn(numNueron(n+1),numNueron(n));
    b{n} = randn(numNueron(n+1),1);
end
%load existing neural network parameters
%load NeuralNetwork9 W b
load NeuralNetwork5_3 W b
%load NeuralNetwork7_3 W b
network.W = W;
network.b = b;
%% Input set
C=[1 0 0 ;  -1 0 0 ;   0 1 0 ;   0 -1 0 ;   0 0 1 ;   0 0 -1];d=[1;1;1;1;1;1];
%C=[1 0 0 0 0;  -1 0 0 0 0;   0 1 0 0 0;   0 -1 0 0 0;   0 0 1 0 0;   0 0 -1 0 0;  0 0 0 1 0;  0 0 0 -1 0; 0 0 0 0 1; 0 0 0 0 -1];d=[1;1;1;1;1;1;1;1;1;1];

%C=[1 0;  -1 0;   0 1 ;   0 -1;];d=[1;1;1;1];

input = Polyhedron(C,d);
%% Compute reachSet of neural network
tic
reachSet = outputNetwork(network,input);
toc
figure;
plot(reachSet)
hold on

figure;
plot(reachSet)
hold on
e=0.05;

for x1 = -1:e:1
    for x2 = -1:e:1
        %for x3 = -1:e:1
        %x=[x1;x2;x3];
        x=[x1;x2];
        for i = 1:1:length(numNueron)-1
            if strcmp(network.activefun(i),'poslin') 
                x = poslin(W{i}*x+b{i});
            elseif strcmp(network.activefun(i),'purelin')
                x = purelin(W{i}*x+b{i});
            end
        end
        if length(x) ==3
            plot3(x(1),x(2),x(3),'*')
        elseif length(x) == 2
            plot(x(1),x(2),'*')        
        end
        %end
    end
end

%save NeuralNetwork W b


figure;plot(reachSet,'linestyle','none','color','green')
hold on
C=[1 0;  -1 0;   0 1 ;   0 -1;];d=[1;1;11;-9];
X0 = Polyhedron(C,d);
h0=plot(X0,'linestyle','none','color','green');
C=[1 0;  -1 0;   0 1 ;   0 -1;];d=[1;1;1;1];
X1 = Polyhedron(C,d);
h1=plot(X1,'linestyle','none','color','b');
C=[1 0;  -1 0;   0 1 ;   0 -1;];d=[1;1;6;-4];
X2 = Polyhedron(C,d);
h2=plot(X2,'linestyle','none','color','red');
legend([h0,h1,h2],'Reachable Set','Unsafe Region 1','Unsafe Region 2')





