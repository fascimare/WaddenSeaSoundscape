% Find relevant times

zID2 = zID(zID(:,2)==2);
Idx = zeros(length(zID2),1);
for n = 1:length(zID2)
    Idx(n) = find(MTT(:) == zID2(n));
end

%Isolate cluster 2

labels2 = allcallID(Idx);

%Isolate cluster 3
zID3 = zID(zID(:,2)==3);
Idx3 = zeros(length(zID3),1);
for n = 1:length(zID3)
    Idx3(n) = find(MTT(:) == zID3(n));
end

labels3 = allcallID(Idx3);