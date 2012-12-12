function [edge orient] = jumpedges(ir,thresholds)

if nargin < 2
    thresholds = [0.90 0.99];
end

thresholds = sort(thresholds);

f = [-1 0 1];
gx=imfilter(ir,f,'replicate','same','conv');
gy=imfilter(ir,f','replicate','same','conv');
edge = sqrt(gx.^2 + gy.^2);

orient = atan(gy./gx)*180/pi;%+90;
orient(orient<0) = orient(orient<0)+180;


% edge = edge./ir;
% edge(isnan(edge)) = 0;
% edge(isinf(edge)) = 0;


h=[];%zeros(size(edge,1)*size(edge,2));
for i=1:size(edge,1);
    h=[h (edge(i,1:size(edge,2)))];
end

h=hist(h,1000);
h=h/sum(h);
%figure;plot(h);hold on;



sumh = 0;
thresmin=0;thresmax=size(h,2);
for i=1:size(h,2);
    sumh=sumh+h(i);
    if sumh > thresholds(2)
        thresmax=i;
        break;
    elseif sumh > thresholds(1) && thresmin==0;
        thresmin=i;
    end
    
end
dif = (max(max(edge))-min(min(edge)));
thresmin=thresmin/size(h,2)*dif;
thresmax=thresmax/size(h,2)*dif;

edge = (edge-thresmin)/(thresmax-thresmin);

edge(edge>1)=1;
edge(edge<0)=0;



end