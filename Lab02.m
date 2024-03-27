s = tf('s');

%% 1
H = 1/(s*(s+2)*(s+4));
generatePlots(H)

%% 2
H = -0.1*(1-2*s)/(s*(s+0.2)*(s+1));
generatePlots(H)

%% 3
H = 1/(s*s*(s+3));
generatePlots(H)

%% 4
H = 2*(1+0.5*s)/((1+s)*(1-s)^2);
generatePlots(H)

%% 5
H = (s*s+1)/((s-2)*(s+2)*(s+4));
generatePlots(H)

%% 7
H = (s-2)/((2+s)*(s*s+1));
generatePlots(H)