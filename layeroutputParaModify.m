function outputPY = layeroutputParaModify(layer,input,activefun)
%Input set is described by a polytope ARRAY Cx < d
%if x = c+sum(a_i.v), input.C*a <d, it means c = 0, V = I so that x = a 
%input.d, where input.c is the center, input.V are the basis vectors, input.C and input.d are the linear constraints 
%Single layer is described by weight layer.W and bias layer.b;
%Output is described a polytope ARRAY PY defined on x(n) > 0, for all n and
%surfaces x(i) >0, x(n) =0, i \ne n
%activefun = poslin or purelin;
numPoly = length(input);
numNeuron=length(layer.b);
emptySet = Polyhedron([],[]); % empty set
outputPY = emptySet;
if strcmp(activefun,'purelin')
    parfor i = 1:1:numPoly
        outputPY(i) = layer.W*input(i)+layer.b;
    end
elseif strcmp(activefun,'poslin') 
    indexProj = dec2bin(0:2^numNeuron -1);
    indexProjMatrix = zeros(2^numNeuron,numNeuron);
    parfor i = 1 : 2^numNeuron
        temp = indexProj(i,:);
        indexProjMatrix(i,:) = str2num(temp(:))';
    end
    indexProjMatrix(2^numNeuron,:) = [];
    indexProjMatrix(1,:) = [];
    I = eye(numNeuron);
    parfor j = 1:1:numPoly*numNeuron %Determine the active and inactive neurons, get the activeVector in which 1 means active, 0 means inactive, 2 means both
    i = ceil(j/numNeuron);  %index of input polytope
    n = j-(i-1)*numNeuron;  %index of neuron
    [l,m] = size(input(i).H);
    f = layer.W(n,:);
    A = input(i).H(:,1:m-1);
    b = input(i).H(:,m);
    Ae = input(i).He(:,1:m-1);
    be = input(i).He(:,m);
    %lb = -inf(m-1,1);
    %ub = inf(m-1,1);
    options = optimset('Display','none');
    xmin = linprog(f,A,b,Ae,be,[],[],[],options);
    ymin = f*xmin+layer.b(n)
    if ymin > 0
        activeVector(j) = 1;
    else
        xmax = linprog(-f,A,b,Ae,be,[],[],[],options);
        ymax = f*xmax+layer.b(n)
        if ymax <= 0
            activeVector(j) = 0;
        else
            activeVector(j) = 2;
        end
    end
end

for i = 1:1:numPoly%Convert the vector to matrix and remove the abundant line
    for j = 1:1:numNeuron
        n = (i-1)*numNeuron+j;
        activeMatrix(i,j) = activeVector(n);
    end
    inactiveIndex =  find(activeMatrix(i,:) == 0);
    if length(inactiveIndex) > 0 %There exists neuron always = 0;
        xflag(i) = 0; % xflag = 0 means all x>0 does not hold for ploytope i;
    else
        xflag(i) = 1;
    end
    row2Remove = 0;
    for k = 1:1:length(inactiveIndex)
        row2Remove = [row2Remove;find(indexProjMatrix(:,inactiveIndex(k)) == 1)];
    end
    activeIndex = find(activeMatrix(i,:) == 1);
    if length(activeIndex) > 0 %There exists neuron always = 0;
        x0flag(i) = 0; % x0flag = 0 means all x<0 does not hold for ploytope i;
    else
        x0flag(i) = 1;
    end
    for k = 1:1:length(activeIndex)
        row2Remove = [row2Remove;find(indexProjMatrix(:,activeIndex(k)) == 0)];
    end
    row2Remove(1) = [];
    indexProjMatrixTemp = indexProjMatrix;
    indexProjMatrixTemp(row2Remove,:) = [];
    indexProjMatrixModify{i} = indexProjMatrixTemp;
    [l,m] = size(indexProjMatrixTemp);
    N(i) = l;
    sumN(i) = sum(N); 
end
    %This part is computing the output y = max(0,W*x+b) with respect to region
    %of all x > 0
    parfor i = 1:1:numPoly
        if xflag(i) == 1
            [n,m] = size(input(i).H);
            V = eye(m-1);
            %compute the output of one layer through y = W*x+b 
            c = layer.b;   %new center c 
            V = layer.W*V;     %new basis vector V
            %add new constraint of x > 0 to C*a < d, i.e, add c+V*a > 0 => -V*a < c
            C = [input(i).H(:,1:m-1);-V];  %new C for constraints on a
            d = [input(i).H(:,m);c];   %new d for constraints on a 
            Ce = input(i).He(:,1:m-1);
            de = input(i).He(:,m);
            PX = Polyhedron('A',C,'b',d,'Ae',Ce,'be',de);
            if PX.isEmptySet() == 0
                PY(i) = V*PX+c;
            else
                PY(i) = emptySet;
            end
        else
            PY(i) = emptySet;
            continue;
        end
    end
    %This part is computing the output y = max(0,W*x+b) with all x(i) < 0, i =
    %0 to number of neurons in the layer, num.      
   % N = length(find(activeVector == 2));
    %N = numPoly*(2^numNeuron-2);
    %i=1;
    parfor j = 1:1:sum(N)
       %if j > sum(N(1:i))
       %    i=i+1;
       %end
       for k = 1:1:length(sumN)
           if sumN(k) >= j
               i = k;
               break;
           end
       end
       n = j-sum(N(1:i-1));
        [l,m] = size(input(i).H);
        V = eye(m-1);
        %compute the output of one layer through y = W*x+b 
        c = layer.b;   %new center c 
        V = layer.W*V;     %new basis vector V
        %compute the region of W*x+b with [W*x+b](n) < 0, i.e.,c(n)+V(n,:)*a < 0 =>
            %V(n,:)*a < -c(n)
        zeroIndex = find(indexProjMatrixModify{i}(n,:) == 0);
        C0 = input(i).H(:,1:m-1);
        d0 = input(i).H(:,m);
        Proj = I;
        for s = 1:1:length(zeroIndex)            
            C0 = [C0;V(zeroIndex(s),:)];
            d0 = [d0;-c(zeroIndex(s))];
            Proj(zeroIndex(s),zeroIndex(s)) = 0;
        end
        c0 = Proj*c;
        V0 = Proj*V;
        C0 = [C0;-V0];
        d0 = [d0;c0];
        Ce0 = input(i).He(:,1:m-1);
        de0 = input(i).He(:,m);
        PX0 = Polyhedron('A',C0,'b',d0,'Ae',Ce0,'be',de0);
        if PX0.isEmptySet() == 0
           PY0(j) = V0*PX0+c0;
        else
           PY0(j) = emptySet;
        end
    end
     
    outputPY = [PY,PY0];
    L = length(outputPY);
    i = 1;
    emptyIndex = 0;
    parfor i = 1:1:L
         if outputPY(i).isEmptySet() == 1;
            emptyIndex = [emptyIndex,i];
         end
    end
    emptyIndex(1) =[];
    outputPY(emptyIndex) = [];
else
    outputPY = 'Invalid.';
    fprintf('This function is only valid for poslin and purelin.');
end









