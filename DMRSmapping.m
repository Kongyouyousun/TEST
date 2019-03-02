function [grid,dmrsindices] = DMRSmapping(grid,DMRS,port)
k = [1 2];
if mod(port,2)==1
    wf = [1 1];
else
    wf = [1 -1];
end
if port>=6
    wt = [1 1];
else
    wt = [1 -1];
end
switch port
    case 1
        tri = 0;
    case 2
        tri = 0;
    case 3
        tri = 2;
    case 4
        tri = 2;
    case 5
        tri = 4;
    case 6
        tri = 4;
    case 7
        tri = 0;
    case 8
        tri = 0;
    case 9
        tri = 2;
    case 10
        tri = 2;
    case 11
        tri = 4;
    otherwise
        tri = 4;
end
dmrsindices = [];
lo = 3;% DMRS放在OFDM3,4的位置上
l = [lo lo+1];
for len = 1:length(DMRS)/2
    n = 6*(len-1);
    for i = 1:2
        for j = 1:2
            grid(n+tri+k(i),l(j))=wf(i)*wt(j)*DMRS(2*(len-1)+k(i));
        end 
        dmrsindices = [dmrsindices;n+tri+k(i)+(lo-1)*size(grid,1)];
    end
end
dmrsindices = [dmrsindices;dmrsindices+size(grid,1)];