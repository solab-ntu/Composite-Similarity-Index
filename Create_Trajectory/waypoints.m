function [waypoint,require_velocity] = waypoints(sys_par,case_name)
path = case_name;
switch(path)
    %% =========== Circle ===============
    % Needed system parameter: radius r
    case 'circle'
        r = sys_par(1);
        t = linspace(0,10,10000);
        x = r*cos(2*r*t-pi/2);
        y = r*sin(2*r*t-pi/2);
        x_ini = linspace(0,r,10);
        y_ini = zeros(1,length(x_ini));
        x = [x_ini x+r];
        y = [y_ini y+r];
        waypoint = [x' y'];
        require_velocity = sys_par(2);
%         plot(waypoint(:,1),waypoint(:,2));
        
        %% ========= DLane ===============
        % Need parameters: [point1, point2,0,zoom_rate,zoom_long]
    case 'DLane'
        require_velocity = 1;
        a1 = sys_par(1);
        a2 = sys_par(2);
        b = sys_par(3);
        zoom_rate = sys_par(4);
        zoom_long = sys_par(5);
        %zoom_rate should be <=zoom_long
        a1 = a1/10;
        a2 = (10-a2)/10;
        b = b/10;
        store_data = [];
        t = linspace(0,1,99);
        pt1 = [0 0;a1 b;a2 (1-b);1 1]*5;
        pt2 = [0 1;1-a2 1-b;1-a1 b;1 0]*5;
        for ii = 1:2
            eval(['pt=pt',num2str(ii),';']);
            pts = kron((1-t).^3,pt(1,:)') + kron(3*(1-t).^2.*t,pt(2,:)') + kron(3*(1-t).*t.^2,pt(3,:)') + kron(t.^3,pt(4,:)');
            pts(1:2,end+1)=pts(1:2,end)+pts(1:2,end)-pts(1:2,end-1);
            ptsX(ii,:) = pts(1,1:length(pts));
            ptsY(ii,:) = pts(2,1:length(pts));
            store_data(2*ii-1,:) = pts(1,:);
            store_data(2*ii,:) = pts(2,:);
        end
        % end
        points = [store_data(1,:) linspace(store_data(1,end),store_data(1,end)+store_data(3,1)+2,20)  store_data(3,:)+store_data(1,end)+2;
            store_data(2,:) ones(1,20)*store_data(2,end) store_data(4,:)];
        over_dis = 3;
        x_points = [linspace(0,over_dis,50) points(1,:)+over_dis linspace(points(1,end)+over_dis,points(1,end)+10*over_dis,500)]*zoom_long;
        y_points = [zeros(1,50) points(2,:) ones(1,500)*points(2,end)]*zoom_rate;
        waypoint = [x_points' y_points'];
        %% =========== Chirp Sine ==================
    case 'chirp'
        require_velocity = 1;
        % w = [0.05,0.1,0.15,0.2,0.25]
        % A = [0.5,0.6,0.7,0.8,0.9,1]
        % gr = 0.01;
        % Needed Parameters: angular frequency omega,amplitude A,growing
        % rate gr
        w = sys_par(1);
        A = sys_par(2);
        gr = sys_par(3);
        t = linspace(0,500,10000);
        y = A*sin((0.5+0.05*exp(gr.*t))*w.*t);
        waypoint = [(t*0.2)' y'];
    %%  ==========Inverse Chirp Sine ================
        case 'inv_chirp'
            % w = [0.04,0.08,0.12,0.16,0.2]
            % A = [0.5,0.6,0.7,0.8,0.9,1]
            % gr = 0.1
        require_velocity = 1;
        sim_time = 30;
        % Needed Parameters: angular frequency omega,amplitude A,growing
        % rate gr
        w = sys_par(1);
        A = sys_par(2);
        gr = sys_par(3);
        t = linspace(0,150,1500);
        y = A*sin((1-0.00000002*exp(gr.*t))*w.*t);
        t = (t*0.2)'+7;
        tk = linspace(0,7,100);
        t = [tk';t];
        y = [zeros(1,100) y];
        
        waypoint = [t y'];
    %% ==========Diverge Sine =====================
    
  
  %% ==========Fixed-Steered Heading ===============
    case 'fix_steer'
  heading = sys_par(1);% heading angle wrt global coordinate
  require_velocity = sys_par(2);
  x1 = linspace(0,7,100);
  y1 = zeros(1,100);
  x2 = linspace(0,50,100);
  y2 = tan(heading*pi/180)*x2;
  x = [x1 x2+7];
  y = [y1 y2];
  waypoint = [x' y'];
  
  %% =========8 ==============================
  
  
  
  
  
end

% plot(waypoint(:,1),waypoint(:,2));
% axis([0,50,0,50]);
% hold on