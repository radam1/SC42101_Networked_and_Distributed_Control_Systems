%% Assignment 1 for Networked and Distributed Control
clear; close all 
%% Part 0: Set up constants/system 
%w/ studnent # 6149022 
a = 6; 
b = 4; 
c = 2; 

A = [0, 0.5-c; 0.2 + a - b, -1]
B = [1;0]

%% Question 1: Design Linear Continuous-Time Controller
sys_c = ss(A, B, eye(2), [0;0]);

%Check step response for ranges of h/tau to check
si = stepinfo(sys_c); 
si.SettlingTime %8.7s and 7.6s respectively

% Controlling using the derived K for pole placement
K = [3, -0.5909]; 

%Controlled System A matrix
A_c = A - B * K; 
eig(A_c) %Double check if you have the right eigenvalues

% Construct the EXACT discrete-time System(Find F(h), G(h)) using symbolic
% math
syms h
A_tilde = [A, B; 0,0,0];
%USE EXPM INSTEAD OF EXP!!!!!
eat = expm(A_tilde .* h)

%Get symbolic expressions for F & G
F = eat(1:2, 1:2)
G = eat(1:2, 3)

%double check by plugging in h=1
F_test = subs(F, h, 1)
G_test = subs(G, h, 1)
discrete_test = c2d(sys_c, 1)

%Check range of sampling times where system is stable:
max_h = 10; 
sample_range = 0.01:0.01:max_h;
system_eigens = zeros(1, length(sample_range)); 
for i = 1:length(sample_range)
    %Get specific sample time
    h = sample_range(i); 
    %Discretize system at the given sampling time
    dsys_i = c2d(sys_c, h); 
    F = dsys_i.A; G = dsys_i.B; 
    dsys_cl = F - G * K; 
    system_eigens(i) = max(abs(eig(dsys_cl))); 
end

f = figure(); 
%Plot gathered sampling times
plot(sample_range, system_eigens, "b")
hold on 
plot([0, max_h], [1,1], "r:")
title("Spectral Radius over Various Sampling Times", FontSize=18)
xlabel("Sampling Time [s]", FontSize=16)
ylabel("Spectral Radius [-]", FontSize=16)
legend("Spectral Radius", "Stability Threshold", FontSize=16)

exportgraphics(f, "Figures/q1_spectral_radius.pdf") 

%Save the workspace(For debugging)
%save("a1")
%% Question 2: Derive stability for tau between 0-h

%Have derived the functions in report Appendix A. 
function [F, G] = get_delayed_dtsys(h, tau)
f = 1.7464; 
Fx = exp(-0.5 * h) * [cos(f*h) + 0.2863 * sin(f*h), -0.859*sin(f*h); 1.26*sin(f*h), cos(f*h) - 0.2863 * sin(f*h)];

Fu = [exp(-0.5 * h) * (0.486*sin(f * h) - 0.303 * cos(f * h)) - exp(-0.5 * (h-tau)) * (0.486*sin(f * (h-tau)) - 0.303 * cos(f * (h-tau)));...
    exp(-0.5 * h) * (-0.191*sin(f * h) - 0.667 * cos(f * h)) - exp(-0.5 * (h-tau)) * (-0.191*sin(f * (h-tau)) - 0.667 * cos(f * (h-tau)))];

G1 = [exp(-0.5 * (h-tau)) * (0.486*sin(f * (h-tau)) - 0.303 * cos(f * (h-tau))) + 0.303;...
    exp(-0.5 * (h-tau)) * (-0.191*sin(f * (h-tau)) - 0.667 * cos(f * (h-tau))) + 0.667];

F = [Fx, Fu; 0, 0, 0]; 

G = [G1; 1];
end

%choose suitable k_u-1 for controller
K_1b = [K,0];
hs = 0.01:0.01:max_h;
spec_radius = zeros(size(hs, 2)); 

%iterate through taus and hs and plot
for h_i=1:length(hs)
    h = hs(h_i); 
    %to make sure tau never goes above h
    h_i_range = length(0.01:0.01:h); 
    for tau_i=1:h_i_range
        tau = hs(tau_i); 
        [Fi, Gi] = get_delayed_dtsys(h, tau); 
        d_sys_cl2 = Fi - Gi * K_1b; 
        spectral_radius = max(abs(eig(d_sys_cl2))); 
        spec_radius(tau_i, h_i) = spectral_radius; 
    end
end

%Plotting Mesh w/ ISO view
f = figure(); 
[X,Y] = meshgrid(hs, hs);
mesh(X, Y, spec_radius)
xlabel("h [s]", FontSize=16)
ylabel("\tau [s]", FontSize=16)
zlabel("Spectral Radius", FontSize=16)
title("Surface Plot of Closed-Loop Spectral Radius", FontSize=18)
view([-45, 45])
exportgraphics(f, "Figures/q2_mesh_iso_specradius.pdf")

%Plotting Surface
f = figure(); 
[X,Y] = meshgrid(hs, hs);
Z = spec_radius;

% Clip values above 1
Z_clipped = min(Z, 1); 

% Create base colormap and add black at the end
cmap = jet(256); 
cmap(end+1,:) = [0 0 0]; % add black

% Create index data: map values >1 to a separate color (black)
Z_index = round(rescale(Z_clipped, 1, 256)); 
Z_index(Z > 1) = 257; 
Z_index(Z == 0) = 257;

% Plot the surface
surf(X, Y, Z, Z_index, 'EdgeColor', 'none')
colormap(cmap)
c = colorbar(); 

colorTitleHandle = get(c,'Title');
set(colorTitleHandle ,'String','Spectral Radius');

c.Ticks = linspace(30, 257, 6);

labels = string(round(linspace(min(min(Z(Z~=0))), 1, 6), 2)); 
labels = [labels(1:5), ">1"]; 
c.TickLabels = labels;

xlabel("h [s]", FontSize=16)
ylabel("\tau [s]", FontSize=16)

title("Heatmap of Closed-Loop Spectral Radius", FontSize=18)
view(2) % for a top-down heatmap view
exportgraphics(f, "Figures/q2_surf_top_specradius_10.pdf")

%to check in a more restricted area
xlim([0,3])
ylim([0,3])
title("Reduced Heatmap of CL Spectral Radius", FontSize=18)
exportgraphics(f, "Figures/q2_surf_top_specradius_3.pdf")
% Selecting sampling interval h to guarentee stability w/ no delay τ
h_22 = 0.4; 
taus = linspace(0, h_22); 
spec_22 = zeros(size(taus));
%Plotting the stability thresholds of 

for i=1:length(taus)
    tau_i = taus(i);
    [Fi, Gi] = get_delayed_dtsys(h_22, tau_i); 
    d_sys_cl2 = Fi - Gi * K_1b; 
    spectral_radius = max(abs(eig(d_sys_cl2)));
    spec_22(i) = spectral_radius; 
end

% Now try with redesigned controller: 
spec_22b = zeros(size(taus));
K_2b = [2, -0.5, 0.25]; 
for i=1:length(taus)
    tau_i = taus(i);
    [Fi, Gi] = get_delayed_dtsys(h_22, tau_i); 
    d_sys_cl2 = Fi - Gi * K_2b; 
    spectral_radius = max(abs(eig(d_sys_cl2)));
    spec_22b(i) = spectral_radius; 
end

%Plotting
f = figure(); 
plot(taus, spec_22)
hold on
plot(taus, spec_22b)
plot([taus(1), taus(end)], [1,1])
hold off
xlabel("\tau [s]", FontSize= 16)
ylabel("Spectral Radius \rho [-]", FontSize= 16)
title(sprintf("Spectral Radius of Two Controllers at h = %0.3fs", h_22), FontSize= 18)
legend("Original Controller", "Redesigned Controller", "Stability Threshold", FontSize= 16)
legend('Location','northwest')
exportgraphics(f, "Figures/q2_controller_comparison.pdf")

% Save the workspace(For debugging)
% save("a1")
%% Question 3
%First redesign initial controller to get poles at -1 and -3
p = [-1, -3]; 

K_3 = place(A, B, p);

%double check to make sure it worked
eig(A - B * K_3)

%Part 1: Deriving the exact discrete-time model. Done in Appendix

%Part 2: Finding ranges of h and tau for stability

%First create the function to get the full closed-loop system 
function A_cl = get_delayed_dtsys_q3(h, tau, K1, K2)
f = 1.7464; 
Fx = exp(-0.5 * h) * [cos(f*h) + 0.2863 * sin(f*h), -0.859*sin(f*h);...
    1.26*sin(f*h), cos(f*h) - 0.2863 * sin(f*h)];

Fu = [exp(-0.5 * h) * (0.486*sin(f * h) - 0.303 * cos(f * h)) - exp(-0.5 * (h-tau)) * (0.486*sin(f * (h-tau)) - 0.303 * cos(f * (h-tau)));...
    exp(-0.5 * h) * (-0.191*sin(f * h) - 0.667 * cos(f * h)) - exp(-0.5 * (h-tau)) * (-0.191*sin(f * (h-tau)) - 0.667 * cos(f * (h-tau)))];

G1 = [exp(-0.5 * (h-tau)) * (0.486*sin(f * (h-tau)) - 0.303 * cos(f * (h-tau))) + 0.303;...
    exp(-0.5 * (h-tau)) * (-0.191*sin(f * (h-tau)) - 0.667 * cos(f * (h-tau))) + 0.667];

F = [Fx, Fu; 0, 0, 0]; 

G = [G1; 1];

A_1 = Fx - G1 * K1;
A_2 = -Fu * K1 - G1 * K2;
A_3 = -Fu * K2;
A_cl = [A_1, A_2, A_3; ... 
       eye(4), zeros(4,2)]; 
end

% At first, have both of the gain matrices be the derived matrices that put
% the poles at -1 and -3
K1_init = K; K2_init = K_3; 

%Using these gain matrices, check stability over a range of h and tau 
hs = 0.01:0.01:max_h;
spec_radius_q3 = zeros(size(hs, 2)); 
%iterate through taus and hs and plot
for h_i=1:length(hs)
    h = hs(h_i); 
    %to make sure tau never goes above h
    h_i_range = length(0.01:0.01:h); 
    for tau_i=1:h_i_range
        tau = hs(tau_i); 
        A_cl = get_delayed_dtsys_q3(h, tau, K1_init, K2_init); 
        spectral_radius = max(abs(eig(A_cl))); 
        spec_radius_q3(tau_i, h_i) = spectral_radius; 
    end
end

f = figure(); 
[X,Y] = meshgrid(hs, hs);
Z = spec_radius_q3;

% Clip values above 1
Z_clipped = min(Z, 1); 

% Create base colormap and add black at the end
cmap = jet(256); 
cmap(end+1,:) = [0 0 0]; % add black

% map values >1 to black
Z_index = round(rescale(Z_clipped, 1, 256)); 
Z_index(Z > 1) = 257; 
Z_index(Z == 0) = 257;

% Plot the surface
surf(X, Y, Z, Z_index, 'EdgeColor', 'none')

%set colormap/colorbar
colormap(cmap)
c = colorbar(); 
colorTitleHandle = get(c,'Title');
set(colorTitleHandle ,'String','Spectral Radius');
c.Ticks = linspace(70, 257, 6);
labels = string(round(linspace(min(min(Z(Z~=0))), 1, 6), 2)); 
labels = [labels(1:5), ">1"]; 
c.TickLabels = labels;

%Rows(y's) are taus, collumns(x's) are h's
xlabel("h [s]", FontSize= 16)
ylabel("\tau [s]", FontSize= 16)

title("Heatmap of Closed-Loop Spectral Radius", FontSize= 18)
view(2) % for a top-down heatmap view

exportgraphics(f, "Figures/q3_surf_top_specradius_w_delays.pdf")
% Now redesign a controller to extend the stability robustness for tau =
% h/2 

K1_q3_p3 = K1_init;
K1s = -5:0.1:5; 
K2s = -5:0.1:5; 
hs = 0.01:0.01:max_h;
K2_best = [0, 0];
best_range = 0; 
for idx_1 =1:length(K1s)
    for idx_2 = 1:length(K2s) 
        %Writing code to maximize the range of 
        K2_q3_p3 = [K1s(idx_1), K2s(idx_2)];
        spec_radius_q3_p3 = zeros(1, length(hs)); 
        %iterate through taus and hs and plot
        for h_i=1:length(hs)
            h = hs(h_i); 
            %to make sure tau never goes above h
            h_i_range = length(0.01:0.01:h); 
            tau = h/2; 
            A_cl = get_delayed_dtsys_q3(h, tau, K1_q3_p3, K2_q3_p3); 
            spectral_radius = max(abs(eig(A_cl))); 
            spec_radius_q3_p3(h_i) = spectral_radius; 
        end

        stability_range = sum(spec_radius_q3_p3 < 1); %# of points in stability region
        if stability_range > best_range
            K2_best = K2_q3_p3; 
            best_range = stability_range; 
        end 
    end
end 

%Now, with the best controller, plot the range of hs that are stable
K2_q3_p3 = K2_best; 
spec_radius_q3_p3 = zeros(1, length(hs)); 
%iterate through taus and hs and plot
for h_i=1:length(hs)
    h = hs(h_i); 
    %to make sure tau never goes above h
    h_i_range = length(0.01:0.01:h); 
    tau = h/2; 
    A_cl = get_delayed_dtsys_q3(h, tau, K1_q3_p3, K2_q3_p3); 
    spectral_radius = max(abs(eig(A_cl))); 
    spec_radius_q3_p3(h_i) = spectral_radius; 
end


f = figure(); 
plot(hs, spec_radius_q3_p3)
hold on 
plot([hs(1), hs(end)], [1,1]) 
hold off
xlabel("h [s]", FontSize= 16)
ylabel("Spectral Radius \rho [-]", FontSize= 16)
legend("CL System Spectral Radius", "Stability Threshold", FontSize= 16)
title("Range of Stable Sampling Times for Controlled System for \tau=\frac{h}{2}", FontSize= 18)
exportgraphics(f, "Figures/q3_controlled_stable_hranges.pdf")

% Save the workspace(For Debugging)
% save("a1")
%% Question 4: Setting LMIs

%Parts 1 and 2 are answered in the notes I have in the ipad

%Part 3: Finding an h1 and h2 to ensure stability
%no time delays, so can just use the F(h) and G(h) from part 2

%For Case 1: Need to compute F and G for both, then use K1 and K2 depending
%on the time (h1h2)^ω
h_range = 0.01:0.01:max_h; 
spectral_radius_case1 = zeros(length(h_range), length(h_range)); 

%Iterate through h1 and h2 to check spectral radius
for h1_i =1:length(h_range)
    h1 = h_range(h1_i); 
    for h2_i = 1:length(h_range)
        h2 = h_range(h2_i);
        %Find the F1 and G1(tau = 0)
        dsys_1 = c2d(sys_c, h1); 
        F1 = dsys_1.A; G1 = dsys_1.B; 
        %Get A1(using K1 because h=h1)
        A1 = (F1 - G1 * K1_q3_p3); 
        
        %Find the F2 and G2(tau = 0)
        dsys_2 = c2d(sys_c, h2); 
        F2 = dsys_2.A; G2 = dsys_2.B;
        %Get A2(using K2 because h=h2)
        A2 = (F2 - G2 * K2_q3_p3); 
        
        %combine to total system and log spectral radius
        A_sys = A1 * A2; 
        spectral_radius_case1(h1_i, h2_i) = max(abs(eig(A_sys))); 

    end
end

f = figure(); 
[X,Y] = meshgrid(h_range, h_range);
Z = spectral_radius_case1;
%surf(X, Y, Z, 'EdgeColor', 'none')

% Clip values above 1
Z_clipped = min(Z, 1); 

% Create a colormap (e.g., jet) and add black at the end
cmap = jet(256); 
cmap(end+1,:) = [0 0 0]; % add black

% Create index data: map values >1 to a separate color (black)
Z_index = round(rescale(Z_clipped, 1, 256)); 
Z_index(Z > 1) = 257; 
Z_index(Z == 0) = 257;

% Plot the surface
surf(X, Y, Z, Z_index, 'EdgeColor', 'none')
colormap(cmap)
c = colorbar(); 
colorTitleHandle = get(c,'Title');
set(colorTitleHandle ,'String','Spectral Radius');

c.Ticks = linspace(2, 257, 6);
label_nums = round(linspace(min(min(spectral_radius_case1)), 1, 6), 2); 

labels = [string(label_nums(1:5)), ">1"];
c.TickLabels = labels;

xlabel("h1 [s]", FontSize= 16)
ylabel("h2 [s]", FontSize= 16)

title("Closed-Loop \rho of (h_1h_2)^{\omega} system wrt h_1 and h_2", FontSize= 18)
view(2) % for a top-down heatmap view

exportgraphics(f, "Figures/q4_case1.pdf")
%For Case 2: Simple Case, only (h1)^ω. This is the same as the previous
%q2. Use controller from last question

%Check range of sampling times where system is stable:
h_range = 0.01:0.01:max_h;
system_eigens = zeros(1, length(h_range)); 
for i = 1:length(h_range)
    %Get specific sample time
    h = h_range(i); 
    %Discretize system at the given sampling time
    dsys_i = c2d(sys_c, h); 
    F = dsys_i.A; G = dsys_i.B; 
    dsys_cl = F - G * K; 
    system_eigens(i) = max(abs(eig(dsys_cl))); 
end

f = figure()
plot(h_range, system_eigens)
hold on 
plot([h_range(1), h_range(end)], [1,1])
hold off
legend("Closed-Loop Spectral Radius", "Stability Threshold")

xlabel("h1 [s]", FontSize= 16)
ylabel("Spectral Radius \rho [-]", FontSize= 16)
title("Closed-Loop \rho of (h_1)^{\omega} system with varying h_1", FontSize= 18)
exportgraphics(f, "Figures/q4_case2.pdf")
save("a1_final")