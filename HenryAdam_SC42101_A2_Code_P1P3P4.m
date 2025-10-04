%% Networked and Distributed Assignment 2
%Declaring Initial Variables
clear; close all
a = 6; 
b = 4; 
c = 2; 

A = [0, 0.5-c; 0.2 + a - b, -1]
B = [1;0]

%% Q1: Omega-Autonomota 
% Item 1: Omega-autonaton can be seen in the notes

% Item 2: Equations of Switched System can be seen in notes

% Item 3: LMI Definitions are in Notes

% Item 4: Plot the stability regions depending on h1 and h2
%no delays, so can use c2d() 

%Get the continuous-time system
ct_sys = ss(A, B, eye(size(A)), zeros(size(B))); 

%define the gain matrices from last assignment: 
K1 = [3.0000,   -0.5909]; K2 = [0.4000,    0.6000]; 

%function for the for loop that will get all of the constituent components
%of the closed-loop matrices
function [F1_h1, F1_h2, F2_h2, F0] = get_system_for_stability(h1, h2, ct_sys, K1, K2)
    dsys_1 = c2d(ct_sys, h1); 
    dsys_2 = c2d(ct_sys, h2); 
    F1 = dsys_1.A; G1 = dsys_1.B; 
    F2 = dsys_2.A; G2 = dsys_2.B; 
    
    F1_h1 = (F1 - G1 * K1);  %state 1
    F1_h2 = (F1 - G1 * K2); %state 2
    F2_h2 = (F2 - G2 * K2); %state 3
    F0 = F2; %state 4
end

%now iterate through several h1 and h2 values to find the 

h1s = 0.01:0.01:10; 
max_spec_radii = zeros(length(h1s),1); 
for i = 1:length(h1s)
    h1 = h1s(i);
    h2 = 2 * h1; 
    [f1h1, f1h2, f2h2, f0] = get_system_for_stability(h1, h2, ct_sys, K1, K2); 
    Fcl1 = f1h1 * f2h2; eig1 = eig(Fcl1); 
    Fcl2 = f1h1 * f0; eig2 = eig(Fcl2); 
    Fcl3 = f1h1 * f2h2 * f1h2 * f2h2; eig3 = eig(Fcl3); 
    Fcl4 = f1h1 * f0 * f1h2 * f0; eig4 = eig(Fcl4); 
    Fcl5 = f1h1 * f2h2 * f1h2 * f0; eig5 = eig(Fcl5); 
    Fcl6 = f1h1 * f0 * f1h2 * f2h2; eig6 = eig(Fcl6); 
    Fcl7 = f2h2 * f1h2; eig7 = eig(Fcl7);
    Fcl8 = f0 * f1h2; eig8 = eig(Fcl8);
    Fcl9 = f2h2 * f1h2 * f0; eig9 = eig(Fcl9);
    Fcl10 = f0 * f1h2 * f2h2; eig10 = eig(Fcl10);
    eigs = [max(abs(eig1)), max(abs(eig2)), max(abs(eig3)), max(abs(eig4)), max(abs(eig5)), max(abs(eig6)), max(abs(eig7)), max(abs(eig8)), max(abs(eig9)), max(abs(eig10))]; 
    %make sure the rows(y) are h2s, the collums(x) are h1s
    max_spec_radii(i) = max(eigs); 
    max(eigs);
end

f = figure();
plot(h1s, max_spec_radii)
hold on 
plot([h1s(1), h1s(end)], [1,1])
hold off
ylabel("Maximum Spectral Radius [-]", FontSize=18)
xlabel("h [s]", FontSize=18)
title("Max Spectral Radius vs h for Switching System", FontSize=18)
legend("Max Spec Radius", "Stability Threshold", FontSize=16)
%exportgraphics(f, "Figures/spec_radius_q1.pdf")
%% Q2: Markov Chain
%Item 1: Markov Chain can be seen in the notes

%Item 2: Equation can be seen in the notes

%Item 3: Stability criteria for Mean-Square Stability can be seen in notes

%Item 4: Unfortunately I couldn't download cvx in MATLAB, so see the python
%code for my solution to the LMIs


%% Q3: Jordan Form
% New dynamics: 
a = 6; 
b = 4; 
c = 2; 

A = [0.3 + a - b, 0; 1, 0.5 + c];

B = [1;0]; 

% Item 1: Pole placement and stability analysis
%desired poles(for some reason place() didnt work but I just solved by
%hand(see notes):
K = [8.8, 21.25]; 

%double check eigenvalues just in case
eig(A-B*K) 

%Since this is a different A matrix, I'll have to re-derive the equations

%check my handwritten work(see along with rest of the math for this in my
%notes)
syms h
Fx = expm(A * h)

%Function to be used in the for loop
function [F, G] = get_delayed_dtsys(h, tau)
Fx = [         exp(2.3 * h),                0;...
      5*exp(2.5*h) - 5 * exp(2.3 * h), exp(2.5 * h)]; 

Fu = [(1/2.3)*(exp(2.3 * h) - exp(2.3 * (h-tau))); ...
       2 * (exp(2.5 * h) - exp(2.5 * (h - tau))) + (5/2.3)*(exp(2.3 * (h - tau)) - exp(2.3 * h))]; 

G1 = [0.435 * exp(2.3*(h -tau)) - 0.4348;...
     2 * exp(2.5 * (h - tau)) - 2.174 * exp(2.3 * (h - tau)) + 0.174];

%Extend the state to include the previous input: 
F = [Fx, Fu; ...
    0,0,0]; 

G = [G1;... 
     1]; 

end


aug_K = [K, 0]; %augment K to act on previous state(but set to zero)
h_range = 0.01:0.01:0.5; 
spec_radii_q2 = zeros(size(h_range, 2), size(h_range, 1)); 

for i = 1:length(h_range)
    h = h_range(i); 
    tau_range = 0.01:0.01:h; 
    for j = 1:length(tau_range)
        tau = tau_range(j); 
        [F, G] = get_delayed_dtsys(h, tau); 
        A_cl = (F - G * aug_K); 
        %make sure h is x, tau is y
        spec_radii_q2(j, i) = max(abs(eig(A_cl))); 
    end
end 

%Create heat plot where all unstable data is black 
[X, Y] = meshgrid(h_range, tau_range); 
Z = spec_radii_q2; 

f = figure();
% Clip values above 1
Z_clipped = min(Z, 1); 

% Create a colormap and add black at the end
cmap = jet(256); 
cmap(end+1,:) = [0 0 0]; 

% map values >1 to a separate color (black)
Z_index = round(rescale(Z_clipped, 1, 256)); 
Z_index(Z > 1) = 257; 
Z_index(Z == 0) = 257;

% Plot the surface
surf(X, Y, Z, Z_index, 'EdgeColor', 'none')
colormap(cmap)
c = colorbar(); 

colorTitleHandle = get(c,'Title');
set(colorTitleHandle ,'String','Spectral Radius');

c.Ticks = linspace(183, 257, 6);

labels = string(round(linspace(min(min(Z(Z~=0))), 1, 6), 2)); 
labels = [labels(1:5), ">1"]; 
c.TickLabels = labels;

xlabel("h [s]", FontSize=16)
ylabel("\tau [s]", FontSize=16)

title("Heatmap of Closed-Loop Spectral Radius", FontSize=18)
view(2) % for a top-down heatmap view
%exportgraphics(f, "Figures/q3_heatmap.pdf")
% Item 2: Jordan Form 

% Item 3: Solving LMIs(See Python Script) 

% Item 4: Comparison of first and second options for graphing. Should find
% that the polytopic overapproximation is a more conservative estimate of
% the stability since it overestimates the matrix space. That means that
% the polytopic estimation should show LESS stable regions than the item 1
% conclusion. 

%% Q4: Event-Triggered Control 
% don't need LMIs, so I can use MATLAB :) 

% Go back to old A,B,etc
a = 6; 
b = 4; 
c = 2; 

A = [0, 0.5-c; 0.2 + a - b, -1];
B = [1;0];

% K from last controller
K = [3.0000,   -0.5909]; 

% Assume Q=I because its positive definite, then use lyap() to solve for
% P(thanks matlab)
Q = eye(2); 
P = lyap((A - B*K)', Q);

t_end = 5; %

% Define the ETC
function [value, isterminal, direction] = ETC_event(t, xt, xk, A, B, K, P, Q, sigma)
    x = [xt; xk];
    middle = [A' * P + P * A + sigma * Q, -P * B * K;...
              -(B * K)' * P, zeros(2)];

    value = x' * middle * x; %when this crosses zero, stop 
    isterminal = 1; % tells ode45 to stop 
    direction = 1; % function should stay below zero, so will flag all crossings where the function is increasing 

end 

% Now do all the actual simulation
init_conditions = [1, 0; 0, 1; 2, 1; -10, -5]; %not exactly sure what to put here so there are just some random ones
sigmas = [0.1, 0.3, 0.5, 0.7, 0.9]; %between zero and one

%arrays to save the means and standard deviations
avg_sampling_times = zeros(1, length(sigmas)); 
std_sampling_times = zeros(1, length(sigmas)); 

for i = 1:length(sigmas)
    sigma = sigmas(i); 
    means = zeros(1, length(init_conditions)); 
    stds = zeros(1, length(init_conditions)); 
    fi = figure(i);
    for ic = 1:length(init_conditions)
        %beginning of the run. make all the stuff you need to do
        %post-processing
        x0 = init_conditions(ic, :)'; 

        %specify the starting conditions of the simulation: 
        t = 0; x = x0; sample_times = t; x_ks = x0; 

        while t(end)<t_end
            x_k = x(:, end); 
            options = odeset('Events', @(t,x) ETC_event(t, x, x_k, A, B, K, P, Q, sigma));
            t_span = [t(end), t_end]; 
            %need to define function in call bcs we need x_k to update
            [t_segment, x_segment, te, xe] = ode45(@(t,x) A*x - B*K*x_k, t_span, x(:,end), options);
            t = [t; t_segment]; %update the current time after each event
            x = [x, x_segment']; % update the current x 
            sample_times = [sample_times, te]; %update the sample times
            x_ki = [ones(1, length(t_segment)) * x_k(1); ones(1, length(t_segment)) * x_k(2)]; 
            x_ks = [x_ks, x_ki]; 
        end
        subplot(2, length(init_conditions), ic)
        plot(t, x(1, :))
        xlabel("time [s]")
        ylabel("x1 [-]")
        title(sprintf("sigma = %0.2f, IC=(%0.2f, %0.2f)", [sigma, x0(1), x0(2)]))
        hold on 
        plot(t, x_ks(1, :))
        hold off 
        legend("CT System", "Samples for Controller")

        subplot(2, length(init_conditions), ic+length(init_conditions))
        plot(t, x(2, :))
        xlabel("time [s]")
        ylabel("x2 [-]")
        hold on 
        plot(t, x_ks(2, :))
        hold off 
        legend("CT System", "Samples for Controller")

        %end of the run. Log mean sampling time for run. 
        sample_periods = diff(sample_times);
        means(ic) = mean(sample_periods); 
        stds(ic) = std(sample_periods); 
        
    end
    %exportgraphics(fi, sprintf("Figures/q4_plots_sigma0_%i.pdf", 100 * sigma))
    avg_sampling_times(i) = mean(means); %avg sampling time across ALL initial conditions
    std_sampling_times(i) = mean(stds); %average standard deviation across all ICs 

end

sigma = ["0.1", "0.3","0.5","0.7","0.9"]';
mean_h = avg_sampling_times'; std_h = std_sampling_times'; 
test_table = table(sigma, mean_h, std_h) 

% Item 3: Testing if system is stable for h=h*
h_star = max(avg_sampling_times); %largest avg time between samples

%this is just a normal sampled-data system w/ constant h, so you can just
%use the normal d2h method

ct_sys = ss(A, B, eye(2), zeros(size(B))); 
d_sys = c2d(ct_sys, h_star); 
F = d_sys.A; G = d_sys.B; 
F_cl = (F - G * K); 
max(abs(eig(F_cl)))
