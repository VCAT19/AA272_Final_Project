% Tyler Weiss AA272 Final Project 3/15/21

% This Script loads the GNSS and POZYX sample data, solves the NR Position
% Alogrithm on the GNSS data (w/ Ionosphere-Free Psuedoranges
% implemented)and then combines the GNSS and POZYX data to simulate a
% Vechiles entering a building.
% It saves the data into a .mat at the end. This .dat will be loaded into a
% different script that will implement an Extended Kalman Filter to
% estimate the position of the quad as it enters the GPS deprived indoor
% state.

clc; clear all; close all;

% Load_Data Sets
Outdoor_Data = readmatrix('gnss_log.csv'); % Load Data into an Array
Outdoor_Data(9704,:) = []; 
Outdoor_Data(9704,:) = [];
Outdoor_Data(9704,:) = [];
Outdoor_Data(9704,:) = [];
Outdoor_Data(9704,:) = [];
Indoor_data = load('h1_front.mat', '-mat');

%Calculate Measured Psuedo Range.
for i = Outdoor_Data(:,1) % Loop over column 1 in Outdoor_Data
    N_w = floor(-(Outdoor_Data(:,4)/604800e9)); % Current GPS Week Number
    t_Rx_hw = Outdoor_Data(:,3)+Outdoor_Data(:,6);% Corrected Reciever Time
    b_hw = Outdoor_Data(:,4)+ Outdoor_Data(:,19); % Adjusted bias in hard clock
    t_Rx_GPS = t_Rx_hw - b_hw; % Reciever time in GPS frame
    t_Rx_w = t_Rx_GPS - (N_w*604800e9); % Reciever time in current week frame
    p_ns = t_Rx_w - Outdoor_Data(:,8); % PsuedoRange in nanoseconds
    pm = p_ns*(299792458/1e9); % PsuedoRange in Meters
    Outdoor_Data_with_PsuedoRange = [Outdoor_Data p_ns pm];
end

% Logic to find unique Time instances and build new arrays based off them,
% then proceed to the next unique time set after running the NR algorithm
% on each unique time instance.
Unique_Time_w_position_EKF = []; % Intialize Array for Solutions outside of Loop
uniqueTimes = unique(Outdoor_Data_with_PsuedoRange(:,2)); % Find Unique times in the second column.
for i = 1 : length(uniqueTimes)
  thisNumber = uniqueTimes(i); % Get a unique Time
  rowsWithThisNumber = Outdoor_Data_with_PsuedoRange(:,2) == thisNumber; % Logical vector.
  for j = rowsWithThisNumber(:,1)
     if rowsWithThisNumber(j,1) == 1 % Find all times that match together.
     Unique_Time_w_data = Outdoor_Data_with_PsuedoRange(j,:); % Pull complete Rows from original dataset
    
     % NR implementation
     x0 = [0, 0, 0]; bu0=0; % Intial guesses
     dx = 10000; dy = 10000; dz = 10000; dbu = 10000; % Intialize variables so While condition is met on first iteration
     x0_x = x0(1); x0_y = x0(2); x0_z = x0(3); 
     %Start Newton Raphson Loop.
     while abs(dx) > .00001 && abs(dy) > .00001 && abs(dz) > .00001 && abs(dbu) > .00001
     
        % Geometry Matrix Loop
        G =[]; % Intialize Geometry Array
        for k = 1: length(Unique_Time_w_data(:,1)) % Calculate Geometry matrix for every data entry at this unique time
            [G_k]=Geometry_Matrix((Unique_Time_w_data(k,22)),(Unique_Time_w_data(k,23)),(Unique_Time_w_data(k,24)), x0_x, x0_y, x0_z);
            G = [G; G_k]; % Add latest result to the Array
        end
        
        % Psuedo Range Loop
        p_t = []; %Intialize theoretical Psuedo Range Array
        for m =1:length(Unique_Time_w_data(:,1)) % Calculate psuedo range for every data entry at this unique time
            [p_t_m]=theoretical_psuedo_range((Unique_Time_w_data(m,22)),(Unique_Time_w_data(m,23)),(Unique_Time_w_data(m,24)),x0_x, x0_y, x0_z,bu0,(Unique_Time_w_data(m,25)));
            p_t = [p_t; p_t_m]; % Add latest result to the Array
        end
        p_i = [];
        
        % Check for any repeated measurements
        for L = 1: length(Unique_Time_w_data(:,5))
            thisSvID = Unique_Time_w_data(L,5);
            rowswiththisSvID = Unique_Time_w_data(:,5) == thisSvID;
            
            %If the measurement is repeated then we must calculate Iono-Free PseudoRange
            if sum(rowswiththisSvID) > 1
                for E = 1:length(rowswiththisSvID)
                    if rowswiththisSvID(E) == 1 
                       % Set frequency and current psuedo range
                       F1 = Unique_Time_w_data(E,16);
                       p1 = Unique_Time_w_data(E,27);
                       % Delete this entry from relevant array's so we can find the other frequency for this SvID in this time instance and pull the data we need
                       rowswithSVID_m = rowswiththisSvID;
                       rowswithSVID_m(E) = [];
                       Unique_Time_w_data_m = Unique_Time_w_data;
                       Unique_Time_w_data_m(E,:) = [];
                       for Q = 1:length(rowswithSVID_m)
                           if rowswithSVID_m(Q) == 1
                               F2 = Unique_Time_w_data_m(Q,16);
                               p2 = Unique_Time_w_data_m(Q,27);
                               break
                           end
                       end
                       pi_i = (((F1^2)/(F1^2-F2^2))*p1) - (((F2^2)/(F1^2-F2^2))*p2);
                       p_i = [p_i; pi_i];
                       break
                    end
                end
            else
                if E > length(rowswiththisSvID) %Handle the case where our first SVID for a time instance is not a repeat
                    E=1;
                end
                %If the measurement is repeated then add the 
                p_i = [p_i; Unique_Time_w_data(L,27)];
            end
        end
        
        %Compute differences between the theoretical and measured psuedoranges. 
        delta_p = p_i - p_t;
     
        %Compute Clock and position Update
        G_T = transpose(G);
        Update_Array = inv((G_T*G))*G_T*delta_p;
        dx = Update_Array(1);
        dy = Update_Array(2);
        dz = Update_Array(3);
        dbu = Update_Array(4);
     
        %Update all states
        x0_x = x0_x + dx; 
        x0_y = x0_y + dy;
        x0_z = x0_z + dz;
        bu0 = bu0 + dbu;
     end
     
     %Obtain some values needed for EKF design
     %Calculate new Tracjectory Time and add it to our array [ms]
     Trajectory_Time = Unique_Time_w_data(1,2) - 1490648213386; % Can hardcode index of 1, since this field will be the same for all rows in the Unqiue_time_with_data array.
     C_N0_std = std(Unique_Time_w_data(:,10)); % [dbHz]
     SvTime_Uncertainty_std = std(Unique_Time_w_data(:,9)); % [ns]
     
     % Append Position to Unique Time array. Define the first position as
     % (0,0,0) 
     Unique_Time_w_position_EKF = [Unique_Time_w_position_EKF; (x0_x+2703950.145499087) (x0_y+4263087.17552874) (x0_z-3885060.0280664838) bu0 Trajectory_Time C_N0_std SvTime_Uncertainty_std];
     
     end
  end
end

Outdoor_End = length(Unique_Time_w_position_EKF);
Indoor_End = length(Indoor_data.pos);

%Get average delta for time in Outdoor_data, then use this as the indoor
%time Delta. Trying this for the indoor data produced 2/9 negative resluts
%so I am questioning the validity of that data set's elapsed time metrics.
prev_time =0;
Indoor_Time_Deltas=[];
for i =1:Outdoor_End
time_i=Unique_Time_w_position_EKF(i,5);
Delta = time_i - prev_time;
prev_time = time_i; % Set up for next iteration
Indoor_Time_Deltas = [Indoor_Time_Deltas; Delta];
end 
Nominal_Ave_Indoor_Time_Delta = mean(Indoor_Time_Deltas);

Set_Over_lap = 5; % Pick how many overlap points we want. 
Mixed_Data =[];
% Merge GNSS Data with the POZYX data.
Indoor_time_i =Unique_Time_w_position_EKF(Outdoor_End-Set_Over_lap,5); % Start the clock for the indoor data wherever the last outdoor time was.
for i = 1 : Indoor_End
    % Add final values of GNSS trajectory to POZXY data then add to array
    Indoor_x_i = Unique_Time_w_position_EKF(end,1) + Indoor_data.pos(i).x;
    Indoor_y_i = Unique_Time_w_position_EKF(end,2) + Indoor_data.pos(i).y;
    Indoor_z_i = Unique_Time_w_position_EKF(end,3) + Indoor_data.pos(i).x;
    Indoor_time_i = Indoor_time_i + (Nominal_Ave_Indoor_Time_Delta/1.1) - ((std(Indoor_Time_Deltas)*rand(1))/2); % Update Time based off average Indoor Time Delta and some noise (Scale to prevent issues with Time)
    
    if i == 1
        Held_data =[];
        for j = 1:Outdoor_End-Set_Over_lap
            Held_data =[Held_data; Unique_Time_w_position_EKF(j,:)];
        end
        Mixed_Data = Held_data;
    end
    
    %Selected 45 and 1000019 to represent poor C/N0 and SvTime_Uncertainity
    %based off worse values found in the GNSS log
    if i <= Set_Over_lap
        Mixed_Data=[Mixed_Data; Indoor_x_i Indoor_y_i Indoor_z_i 0 Indoor_time_i 45 1000019; Unique_Time_w_position_EKF(Outdoor_End-Set_Over_lap+i,:)];
    else
        Mixed_Data=[Mixed_Data; Indoor_x_i Indoor_y_i Indoor_z_i 0 Indoor_time_i 45 1000019];
    end 
end

save('ProjectData.mat', 'Mixed_Data')

% Functions called above
function [G] = Geometry_Matrix (X_k,Y_k,Z_k,x_u,y_u,z_u)
d = sqrt((X_k-x_u)^2 + (Y_k-y_u)^2 + (Z_k-z_u)^2);
XLOS = (X_k-x_u)/d;
YLOS = (Y_k-y_u)/d;
ZLOS = (Z_k-z_u)/d;
G = [-XLOS -YLOS -ZLOS 1];
end

function [p] = theoretical_psuedo_range(X_k,Y_k,Z_k,x_u,y_u,z_u,bu,B)
p = (sqrt((X_k-x_u)^2 + (Y_k-y_u)^2 + (Z_k-z_u)^2))+bu-B;
end