function Sat3DnodeT_v1_f(W1,Ltot,L45,Tv,usunv,uplanv)
    %function for plotting the satellite node surfaces in 3D coloured according
    %to their temperature. 
    %Created by: Antonio Acosta Iborra. 12 March 2021. Carlos III University of
    %Madrid
    %INPUT ARGUMENTS: 
    %W1 [m], lateral side of nodes 1 and 4
    %Ltot [m], length of the satellite (distance from node 1 to node 4.
    %L45 [m], distance between nodes 4 and 5.
    %Tv [K], column or row vector of nodes temperatures Tv=[T1, T2, ....T9]';
    %usunv [-], column or row unitary vector pointing at the Sun
    %uplanv [-], column or row unitary vector pointing at the orbited planet
    %OUTPUT: 
    %3D surface coloured according to the nodes temperatures. 
    %The distance between nodes 4 and 5 is exagenerated to enhance visualization.
    %The long anthena is placed as a reference between nodes 1 and 6. 
    %The Sun direction is the red arrow and planet direction is the blue
    %arrow.
    
    dx=W1;
    dy=W1;
    dz=Ltot;
    dgap=L45*50;
    colormap(jet);
    
    for i = 1:2

        subplot(1,2,i);
        %%%ploting of the external nodes surfaces

        %Node-1 surface
        x=[-dx dx; -dx dx];
        y=[-dy -dy; dy dy];
        z=[0 0; 0 0];
        w=[Tv(1) Tv(1); Tv(1) Tv(1)];
        surface(x,y,z,w);
        hold on
        shading faceted;

        %Node-4 surface
        x=[-dx dx; -dx dx];
        y=[-dy -dy; dy dy];
        z=[dz dz; dz dz];
        w=[Tv(4) Tv(4); Tv(4) Tv(4)];
        surface(x,y,z,w);

        %Node-5 surface
        x=[-dx dx; -dx dx];
        y=[-dy -dy; dy dy];
        z=[dz+dgap dz+dgap; dz+dgap dz+dgap];
        w=[Tv(5) Tv(5); Tv(5) Tv(5)];
        surface(x,y,z,w);

        %Node-6 surface
        x=[dx dx; dx dx];
        y=[-dy dy; -dy dy];
        z=[0 0; dz dz];
        w=[Tv(6) Tv(6); Tv(6) Tv(6)];
        surface(x,y,z,w);

        %Node-7 surface
        x=[-dx dx; -dx dx];
        y=[dy dy; dy dy];
        z=[0 0; dz dz];
        w=[Tv(7) Tv(7); Tv(7) Tv(7)];
        surface(x,y,z,w);

        %Node-8 surface
        x=[-dx -dx; -dx -dx];
        y=[-dy dy; -dy dy];
        z=[0 0; dz dz];
        w=[Tv(8) Tv(8); Tv(8) Tv(8)];
        surface(x,y,z,w);

        %Node-9 surface
        x=[-dx dx; -dx dx];
        y=[-dy -dy; -dy -dy];
        z=[0 0; dz dz];
        w=[Tv(9) Tv(9); Tv(9) Tv(9)];
        surface(x,y,z,w);
        
        %Colorbar on right plot
        if i == 2
            c = colorbar;
            c.Label.String = 'T [K]';
            c.Label.Interpreter = 'latex';
            c.Label.FontSize = 11;
        end
        caxis([200 500])

        %anthena
        Lant6=4*dx;
        Lant7=2*dx;
        Lant8=2*dx;
        Lant9=2*dx;
        plot3([dx dx+Lant6],[0 0],[0 -Lant6],'k');
        plot3([0 0],[dy dy+Lant7],[0 -Lant7],'k');
        plot3([-dx -dx-Lant8],[0 0],[0 -Lant8],'k');
        plot3([0 0],[-dy -dy-Lant9],[0 -Lant9],'k');

        %Sun and Earth directions
        Lsvec=3*dx;
        Lpvec=2*dx;
        vdirv=[-Lant6; -Lant6; -Lant6]/2;
        plot3(vdirv(1),vdirv(2),vdirv(3),'k+');
        usunv = usunv*Lsvec;
        uplanv = uplanv*Lpvec;
        quiver3(vdirv(1),vdirv(2),vdirv(3),usunv(1),usunv(2),usunv(3),'r');
        quiver3(vdirv(1),vdirv(2),vdirv(3),uplanv(1),uplanv(2),uplanv(3),'b');

        %Axes and view
        set(gca,'visible','off');
        set(gca,'xtick',[]);
        if i == 1
            view(3)
        elseif i == 2
            view(37.5, -30);
        end
        axis equal
        axis tight
        set(gcf,'Position',[276.2000  493.8000  771.8000  268.2000]);
        hold off
    end
end