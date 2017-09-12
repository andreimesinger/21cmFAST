% matlab function to plot slice
%
% mostly written by: Dave Spiegel
% usage:


function cm21_image_slice_withps(infile, outfile, dim1, color_bar_flag, expand_factor, varargin)

if length(varargin) > 0
    z = varargin{1};
    xH = varargin{2};
    aveTb = varargin{3};
end
if length(varargin) > 3
    lower_color_limit = varargin{4};
    upper_color_limit = varargin{5};
    ps_filename = varargin{6};
end

fid1 = fopen(infile,'r','n');
dim2=dim1;
dim3=dim1;
variable = fread(fid1,[dim1,dim2*dim3],'real*4');
variable = reshape(variable,dim1,dim2,dim3);
slice = variable(:,24,:);
squeezed_slice = squeeze(slice);
squeezed_slice = squeezed_slice / (1+str2double(z)/10.0)^0.5;

% load power spectrum
[k, PS, Sig] = textread(ps_filename,'%f %f %f');



fontname = 'Charcoal CY';

Xinches = 15;
Yinches = 7;
clf
dpi     = 95;
vscale = 1;
set(gcf,'Position',[5 52 Xinches*dpi Yinches*dpi]);
set(gcf,'PaperPosition',[0.25 0.25 Xinches Yinches])
h = gcf;

h_s1 = subplot(1,2,1);
x = 1:100;
y = 1:100;
imagesc(x,y,squeezed_slice); set(gca,'ydir','normal')
set(gca,'ydir','normal','xtick',[],'ytick',[])
text(103,60,'1 Gpc','fontsize',16,'interpreter','latex','fontweight','bold','fontname',fontname,'rotation',270)

%colormap
%axis equal
%axis tight
axis image

if length(varargin) > 2
    cvec = [lower_color_limit upper_color_limit];
else
    cvec = caxis;
end
caxis(cvec)

%colormap autumn;
l1 = 100;
l2 = 100;
cm_a = autumn(l1);
%colormap winter;
cm_s = winter(l2);
cm = [cm_s; cm_a(end:-1:1,:)];
%cm(end) = [];
cm(l1,:) = [0 0 0];
cm(l1+1,:) = [0 0 0];
%cm = cm*0;
%cm(1,:) = [1 1 1];
colormap(cm)
set(gca,'linewidth',1)

h_s2 = subplot(1,2,2);
loglog(k,PS,'linewidth',1)
set(gca,'xminortick','on','yminortick','on','ytick',logspace(-6,3,10))
xlabel('$\rm k ~ (Mpc^{-1})$','fontsize',16, 'interpreter','latex')
ylabel('$\rm \langle {\delta}\hspace{-0.07in}T_b\rangle^2 ~ \Delta^2_{21}(k) ~ (mK^2)$','fontsize',16,'interpreter','latex')
%ylabel('$\bar{\delta T_b}^2$ \Delta^2_{21}(k) (mK^2)','fontsize',20,'interpreter','latex')
set(gca,'linewidth',1,'fontsize',14);
h_s2_pos = get(h_s2,'position');
%h_s2_pos
set(h_s2,'position',[0.615 0.13 0.325 0.745])
axis([5e-3 2e0 1e-6 2e3])


subplot(h_s1);

%colormap(cm);                                                                                                                                                                                                            
%[m,n] = size(this_array_sorted);                                                                                                                                                                                         
%this_array_sorted_RGB = zeros([size(this_array_sorted), 3]);                                                                                                                                                             
%for ii = 1:m                                                                                                                                                                                                             
%    for jj = 1:n                                                                                                                                                                                                         
%        color_ndx = floor(1 + (this_array_sorted(ii,jj)-cvec(1))/(cvec(end)-cvec(1)) * (num_colors-1));                                                                                                                  
%        this_array_sorted_RGB(ii,jj,:) = cm(color_ndx,:);                                                                                                                                                                
%    end                                                                                                                                                                                                                  
%    if vals(ii) <= n                                                                                                                                                                                                     
%        this_array_sorted_RGB(ii,vals(ii),:) = [0 0 0];                                                                                                                                                                  
%    end                                                                                                                                                                                                                  
%end                                                                                                                                                                                                                      
%hh = image(linspace(-2.5,2.5,501),y,this_array_sorted_RGB);                                                                                                                                                              
%axis image;                                                                                                                                                                                                              
%set(gca,'ytick',y)                                                                                                                                                                                                       




%if we want to print labels on the plot
if length(varargin) > 0
    text(103,97, strcat('$\bf z=', z, '$'), 'fontsize',16,'color','k', 'FontWeight', 'Bold','Fontname',fontname,'HorizontalAlignment','left',  'interpreter','latex')
    text(101,90, strcat('$\rm \bf \langle x_{HI}\rangle_v =', xH, '$'), 'fontsize',16,'color','k', 'FontWeight', 'Bold','Fontname',fontname,'HorizontalAlignment','left', 'interpreter','latex')
    text(101,83, strcat('$\rm \bf \langle {\delta}\hspace{-0.07in}T_b\rangle_v =', aveTb, '$'), 'fontsize',16,'color','k', 'FontWeight', 'Bold','Fontname',fontname,'HorizontalAlignment','left', 'interpreter','latex')
end

%if color bar is wanted
if color_bar_flag > 0
    orig_axes = gca;
    ax_pos = get(gca,'position');
    c_h = axes;
    set(c_h, 'position', [ax_pos(1) (ax_pos(2)-0.04) ax_pos(3) 0.04])
    colors = linspace(cvec(1),cvec(2),1000);
    x = colors;
    y = 0:1;
    imagesc(x,y,colors)
    colormap(cm)
    set(gca,'YTick',[]);
    set(gca,'XTick',[-50:10:50]);
    text(mean(x)-(0.15*diff(cvec)), -2, '{\delta}T_{b} [(1+z)/10]^{-1/2} (mK)', 'fontsize',12,'color','k', 'fontname',fontname,'fontweight','bold')
    axes(orig_axes);
    printstatement = ['print -dtiff ' outfile]
    printstatement
    eval(printstatement);
else
    printstatement = ['print -depsc ' outfile];
    eval(printstatement);

    clf
    c_h = axes;

    set(c_h, 'position', [0.2 0.1 0.6 0.03])
    colors = linspace(cvec(1),cvec(2),100);
    x = colors;
    y = 0:1;
    imagesc(x,y,colors)
    set(gca,'YTick',[]);

    text(mean(x)-(0.017*diff(cvec)),0.8, 'x_{HI}', 'fontsize',12,'color','k', 'fontname',fontname,'fontweight','bold')

    print -depsc colorbar.eps
end


asdf = 1;

