function reachSet = outputNetwork(network,input)

for i = 1:1:length(network.W)
    layer.W=network.W{i};
    layer.b=network.b{i};
    %tic
    if strcmp(network.activefun(i),'poslin') 
        %input = layeroutput(layer,input,'poslin')
        %input = layeroutputPara(layer,input,'poslin')
        input = layeroutputParaModify(layer,input,'poslin')
        %input = quickConvexHull(input)
    elseif strcmp(network.activefun(i),'purelin')
        %input = layeroutput(layer,input,'purelin')
        %input = layeroutputPara(layer,input,'purelin')
        input = layeroutputParaModify(layer,input,'purelin')
    end
    %toc
end
reachSet = input;