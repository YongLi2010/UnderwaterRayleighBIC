function cmap=ni2019_viridis(n)
%NI2019_VIRIDIS Compact viridis approximation matching 2607.07228v1 plots.
if nargin<1, n=256; end
anchors=[68 1 84;72 40 120;62 74 137;49 104 142;38 130 142; ...
    31 158 137;53 183 121;110 206 88;181 222 43;253 231 37]/255;
x=linspace(0,1,size(anchors,1)); xi=linspace(0,1,n);
cmap=interp1(x,anchors,xi,'pchip');
cmap=min(max(cmap,0),1);
end
