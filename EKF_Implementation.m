% Tyler Weiss AA272 Final Project 3/15/21
clc; clear all; close all;
%%
% Load Data
Struct = load('ProjectData', '-mat');
Data = Struct.Mixed_Data;
Measurement_Data = readmatrix('hw4_prob3');

% We are only concerned with the Data around the transition from indoor to
% outdoor localization. The array is built on a hardcoded index based on
% the "Set_Over_lap" variable set in the "Load_And_Generate_Data_Set.m"
% script.
Transition_Data = [];
 for i = 990:1010
     Transition_Data = [Transition_Data; Data(i,1) Data(i,2) Data(i,3) Data(i,4) Data(i,5) Data(i,6) Data(i,7)];
 end 

% Add acceleratometer data from HW4 to our Data Array for use in the EKF.
ax1 = -.0768; ay1 = -.0753; az1 = -.0081; Transition_Data(1,8) = ax1; Transition_Data(1,9) = ay1; Transition_Data(1,10) = az1;
ax2 = -.0506; ay2 = -.0995; az2 = -.0059; Transition_Data(2,8) = ax2; Transition_Data(2,9) = ay2; Transition_Data(2,10) = az2;
ax3 = -.0209; ay3 = -.0218; az3 = -.0222; Transition_Data(3,8) = ax3; Transition_Data(3,9) = ay3; Transition_Data(3,10) = az3;
ax4 = -.0675; ay4 = -.0363; az4 = 0.0518; Transition_Data(4,8) = ax4; Transition_Data(4,9) = ay4; Transition_Data(4,10) = az4;
ax5 = -0.1094; ay5 = 0.0168; az5 = -.0422; Transition_Data(5,8) = ax5; Transition_Data(5,9) = ay5; Transition_Data(5,10) = az5;
ax6 = -.0926; ay6 = -.0987; az6 = .0307; Transition_Data(6,8) = ax6; Transition_Data(6,9) = ay6; Transition_Data(6,10) = az6;
ax7 = -.0611; ay7 = -.0270; az7 = -.00643; Transition_Data(7,8) = ax7; Transition_Data(7,9) = ay7; Transition_Data(7,10) = az7;
ax9 = -.0197; ay9 = -.0007729; az9 = .0781; Transition_Data(9,8) = ax9; Transition_Data(9,9) = ay9; Transition_Data(9,10) = az9;
ax11 = .0677; ay11 = -.0092; az11 = -.00448; Transition_Data(11,8) = ax11; Transition_Data(11,9) = ay11; Transition_Data(11,10) = az11;
ax13 = -.0112; ay13 = -.0169; az13 = -.0033; Transition_Data(13,8) = ax13; Transition_Data(13,9) = ay13; Transition_Data(13,10) = az13;
ax15 = .0300; ay15 = .0763; az15 = .00824; Transition_Data(15,8) = ax15; Transition_Data(15,9) = ay15; Transition_Data(15,10) = az15;
ax17 = .478; ay17 = .0813; az17 = -.112; Transition_Data(17,8) = ax17; Transition_Data(17,9) = ay17; Transition_Data(17,10) = az17;

% Create Truth Data Set. Will be specific to your data so parse your data
% for outliers and then remove those indexes below or implement some code to throw out outliers.
% Essentially this is the trajectory we want the EKF state estimation to
% resemble as it should use the covairance matrcies to not trust the
% outlier points. 
Truth_data = Transition_Data;
Truth_data(8,:)=[];
Truth_data(10,:)=[];
Truth_data(11,:)=[];
Truth_data(12,:)=[];
Truth_data(13,:)=[];
figure(1)
plot3(Truth_data(:,1),Truth_data(:,2),Truth_data(:,3))
grid on
title('"Truth" 3D Trajectory')
figure(2)
plot(Truth_data(:,5),Truth_data(:,1));
grid on
title('"Truth" X Trajectory')
figure(3)
plot(Truth_data(:,5),Truth_data(:,2));
grid on
title('"Truth" Y Trajectory')
figure(4)
plot(Truth_data(:,5),Truth_data(:,3));
grid on
title('"Truth" z Trajectory')
%% Compute standard deviations of accelerations. Will need this values in
%Kalman Filter Implementation.
sigma_store_Ax = [];
sigma_store_Ay = [];
sigma_store_Az = [];
T_store = 0;
T =0;
for A = 1: length(Measurement_Data)
        if Measurement_Data(A,5) ~= 0
           ax_i = Measurement_Data(A,5);
           sigma_store_Ax = [sigma_store_Ax; ax_i];
        end
        if Measurement_Data(A,6) ~= 0
           ay_i = Measurement_Data(A,6);
           sigma_store_Ay = [sigma_store_Ay; ay_i];
        end
        if Measurement_Data(A,7) ~= 0
           az_i = Measurement_Data(A,7);
           sigma_store_Az = [sigma_store_Az; az_i];
        end
end
Sigma_Ax = std(sigma_store_Ax);
Sigma_Ay = std(sigma_store_Ay);
Sigma_Az = std(sigma_store_Az);
Sigma_X = std(Transition_Data(:,1));
Sigma_Y = std(Transition_Data(:,2));
Sigma_Z = std(Transition_Data(:,3));
X_norm = Sigma_X^2/(Sigma_X^2+Sigma_Y^2+Sigma_Z^2);
Y_norm = Sigma_Y^2/(Sigma_X^2+Sigma_Y^2+Sigma_Z^2);
Z_norm = Sigma_Z^2/(Sigma_X^2+Sigma_Y^2+Sigma_Z^2);
Sigma_C_N0 = std(Transition_Data(:,6));
%% EKF implementation

% Pick your R and Q matrices from the different choices below.
Rchoice = 3;
Qchoice = 3;

% Set up some intial variables for the intial iteration
x_prev = [-350; 2170; 2181; 0; 0; 0]; 
P = zeros(6); % States are given to use so initial covariance can be 0's.
state_store = [];
T_store =0;
for i= 1:length(Transition_Data)
     
    %Get some variables we will need 
     if i == 1
        Delta_T =0;
     else
        Delta_T = Transition_Data(i,5)- Transition_Data(i-1,5);
     end
     
     ax = Transition_Data(i,8);
     ay = Transition_Data(i,9);
     az = Transition_Data(i,10);
     Norm_C_N0 = Transition_Data(i,6)/Sigma_C_N0;
     % [x,y,z,b,x_dot,ydot zdot]
     % Define State Model
     F = [1 Delta_T 0 0 0 0;
          0 1 Delta_T 0 0 0;
          0 0 1 Delta_T 0 0;
          0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1];
      
      % Symbolic for use with Jacobian function later
        syms x_s y z x_dot y_dot z_dot
        f = [x_s*Delta_T;
           y*Delta_T;
          z*Delta_T;
          x_dot;
          y_dot;
          z_dot];
      
      X_State = [x_s; y; z; x_dot; y_dot; z_dot];
      
      G = [.5*(Delta_T)^2*ax;
           .5*(Delta_T)^2*ay;
           .5*(Delta_T)^2*az;
           Delta_T*ax;
           Delta_T*ay;
           Delta_T*az;];
     
      x = F*x_prev+G;
     % Define Measurement Model
      H = [1 0 0 0 0 0 
           0 1 0 0 0 0 
           0 0 1 0 0 0];
       
      g = [x_s; y; z];
      
      %Very Basic Q using standard Deviations. This a Q from the robot
      % model. 
      if Qchoice == 1
         Q = [ Sigma_Ax^2 0 0 0 0 0;
            0 Sigma_Ay^2 0 0 0 0;
            0 0 Sigma_Az^2 0 0 0;
            0 0 0 Sigma_Ax^2 0 0;
            0 0 0 0 Sigma_Ay^2 0;
            0 0 0 0 0 Sigma_Az^2];
      end
      
      
      %Now a little fancier. 
      if Qchoice == 2
         Q = [ .25*(Sigma_Ax)^2*(Delta_T)^4 .5*(Sigma_Ax)^2*(Delta_T)^3 0 0 0 0 0;
            0 .25*(Sigma_Ay)^2*(Delta_T)^4 .5*(Sigma_Ay)^2*(Delta_T)^3 0 0 0 0;
            0 0 .25*(Sigma_Az)^2*(Delta_T)^4 .5*(Sigma_Az)^2*(Delta_T)^3 0 0 0;
            0 0 0 .25*(Sigma_b)^2*(Delta_T)^4 .5*(Sigma_b)^2*(Delta_T)^3 0 0;
            0 0 0 .5*(Sigma_Ax)^2*(Delta_T)^3 (Sigma_Ax)^2*(Delta_T)^2 0 0;
            0 0 0 0 .5*(Sigma_Ay)^2*(Delta_T)^3 (Sigma_Ay)^2*(Delta_T)^2 0;
            0 0 0 0 0 .5*(Sigma_Az)^2*(Delta_T)^3 (Sigma_Az)^2*(Delta_T)^2];
      end
      
      % What if we don't trust our dynamics at all.
      if Qchoice == 3
        Q = [100000 0 0 0 0 0;
            0 100000 0 0 0 0;
            0 0 100000 0 0 0;
            0 0 0 100000 0 0;
            0 0 0 0 100000 0;
            0 0 0 0 0 100000];
      end
      
      % Just a little trust
      if Qchoice == 4
        Q = [.9 0 0 0 0 0 ;
            0 .9 0 0 0 0 ;
            0 0 .9 0 0 0 ;
            0 0 0 .9 0 0 ;
            0 0 0 0 .9 0 ;
            0 0 0 0 0 .9 ];
      end
      
      if Rchoice == 1
        R = [Transition_Data(i,1)^2 0 0 
             0 Transition_Data(i,2)^2 0 
             0 0 Transition_Data(i,3)^2 ];
      end 
      
      if Rchoice == 2
        R = [X_norm^2 Norm_C_N0^2 0; 
             Norm_C_N0^2 Y_norm^2 Norm_C_N0^2; 
             0 Norm_C_N0^2 Z_norm^2];
      end
      
      if Rchoice == 3
        R = [Transition_Data(i,1) Transition_Data(i,7) 0; 
             Transition_Data(i,7) Transition_Data(i,2) Transition_Data(i,7); 
             0 Transition_Data(i,7) Transition_Data(i,3)];
      end
      
      if Rchoice == 4
        R = [Transition_Data(i,7) 0 0; 
             0 Transition_Data(i,7) 0; 
             0 0 Transition_Data(i,7)];
      end 
      
      T_store =  T_store + Delta_T;
      [X, P] = E_Kalman_filter(x_prev, F, P, Q, R, g, f, X_State);
      state_store = [state_store; x_prev(1) x_prev(2) x_prev(3) x_prev(4) x_prev(5) x_prev(6) T_store];
      x_prev = X; % Set up for next iteration.
end
figure(5)
plot3(state_store(:,1),state_store(:,2),state_store(:,3))
title('Extended Kalman Filter Trajectory')
grid on
figure(6)
plot(state_store(:,7),state_store(:,1));
grid on
title('EKF X Trajectory')
figure(7)
plot(state_store(:,7),state_store(:,2));
grid on
title('EKF Y Trajectory')
figure(8)
plot(state_store(:,7),state_store(:,3));
grid on
title('EKF Z Trajectory')
function [x, P] = E_Kalman_filter(x, A, P, Q, R, g, f, X_state)
   % Predict
   x = A*x;
   
   A = subs(jacobian(f,X_state),X_state,x);% Need Jacobain of Dynamics Model to linearize.
   A = double(A); %convert from sym to double, greatly reduces runtime
   P = A*P*A' + Q;
   
   H = subs(jacobian(g,X_state),X_state,x); % Need Jacobain of Dynamics Model to linearize.
   H = double(H); %convert from sym to double, greatly reduces runtime
   z = H*x;
   
   % We'll need an Identity matrix in upate.
   length_s = length(x);
   I = eye(length_s); 
   
   % Update
   y = z - H*x;
   K = P*H'*inv(R + H*P*H');
   x = x + K*y;
   P = (I - K*H)*A*(I-K*H)'+(K*R*K');
   y = z - H*x;
end 