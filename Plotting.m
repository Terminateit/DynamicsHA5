function Plotting(t,time, x,y,z,xname,yname,zname, label, number, png_name)
    hold on;
    fig = figure(number);
    set(gcf,'color','w');
    plot(time, double(subs(x,t, time)),time, double(subs(y,t, time)),time, double(subs(z,t, time)))
    ylabel(label, 'FontSize',12);
    xlabel('Time', 'FontSize',12);
    legend(xname, yname, zname,'interpreter','latex','FontSize',15)
    frame = getframe(fig);
    im = frame2im(frame);
    [img,map] = rgb2ind(im,256);
    imwrite(img,map,png_name,'png');
    hold off;
    close all;
end

