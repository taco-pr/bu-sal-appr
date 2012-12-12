function [ nx,ny,nz ] = normals(x,y,z,h)

if nargin < 4
    h=[-1 0 1];
end

%h=[-1 0 0; -1 0 1; 0 0 1];
%h2=fspecial('gaussian',21,10);
%h2=[-1 0 0; 0 0 0; 0 0 1];

%h=conv(h,h2(10,:));

xu=imfilter(x,h,'replicate','same','corr');
xv=imfilter(x,h','replicate','same','corr');
yu=imfilter(y,h,'replicate','same','corr');
yv=imfilter(y,h','replicate','same','corr');
zu=imfilter(z,h,'replicate','same','corr');
zv=imfilter(z,h','replicate','same','corr');

% calculate normals
nx=(yu.*zv - zu.*yv);
ny=(-xu.*zv + xv.*zu);
nz=(xu.*yv - yu.*xv);
% nr=double(sqrt(nx.^2+ny.^2+nz.^2));
% nx=nx./nr;
% ny=ny./nr;
% nz=nz./nr;

h1=[-1 0 0;0 0 0; 0 0 1];
h2=[0 0 1;0 0 0; -1 0 0];
%h=[-1 0 0; -1 0 1; 0 0 1];
%h2=fspecial('gaussian',21,10);
%h2=[-1 0 0; 0 0 0; 0 0 1];

%h=conv(h,h2(10,:));

xu=imfilter(x,h1,'replicate','same','corr');
xv=imfilter(x,h2,'replicate','same','corr');
yu=imfilter(y,h1,'replicate','same','corr');
yv=imfilter(y,h2,'replicate','same','corr');
zu=imfilter(z,h1,'replicate','same','corr');
zv=imfilter(z,h2,'replicate','same','corr');

% calculate normals
nx=nx+(yu.*zv - zu.*yv);
ny=ny+(-xu.*zv + xv.*zu);
nz=nz+(xu.*yv - yu.*xv);
nr=double(sqrt(nx.^2+ny.^2+nz.^2));
nx=nx./nr;
ny=ny./nr;
nz=nz./nr;





end

