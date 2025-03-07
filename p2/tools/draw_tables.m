function draw_tables (compatibility, GT, H)
%-------------------------------------------------------
% University of Zaragoza
% Authors:  J. Neira, J. Tardos
%-------------------------------------------------------
global configuration;

if configuration.step_by_step

ncol=256 ;
mapa=hsv2rgb([linspace(2/3, 0, ncol)' 0.9*ones(ncol,1) ones(ncol,1)]) ;
mapa(1,:)=[1 1 1] ;
mapa = flipdim(mapa,1);

m = colormap;
m(1,:)=[1 1 1] ;

figure(configuration.tables); clf; 
subplot(2,2,1);

colormap(m); 
%table = compatibility.ic(:, compatibility.candidates.features);
table = compatibility.ic;
imagesc(table);

title('INDIVIDUAL COMPATIBILITY');

subplot(2,2,2);
colormap(m);

table = make_table(compatibility, GT);
imagesc(table);
%colorbar;
title('GROUND TRUTH');

subplot(2,2,3);
colormap(m);
table = make_table(compatibility, H);
imagesc(table);
title(configuration.name);

pause
end

function table = make_table (compatibility, H),

table = zeros(size(compatibility.ic));
i = find(H);
if nnz(i)
    j = H(i);
    i = sub2ind(size(table), i, j);
    table(i) = 1;
end
%table = table(:, compatibility.candidates.features);
