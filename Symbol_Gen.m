%%
% Symbol_Gen.m                                          Iniitial: 11/16/20
%
% Description: Generate a dataset of noisy IQ samples for a given
% modulation scheme and interference constellation.
%
% NOTES:
%   - L_ indicates Label (for ML analysis)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all


%% General Paramters
PARAMS.VISUALIZE = 1;          % Set to 1 to show received constellations
PARAMS.NORMALIZE_RESULTS = 1;  % Normalize results to power = 1

%% Simulation Parameters
Len_block = 50;    % Number of symbols per block
Num_blocks = 100;     % Number of blocks observed in dataset
p_int = 0.75;        % Probability of interference being present in a block

%SNR = 100*ones(1,Num_blocks);   % Constant for all blocks
SNR_MIN = 20; SNR_MAX = 100;
SNR = (SNR_MAX-SNR_MIN)*rand(1,Num_blocks) + SNR_MIN;

SIR = 10;                       % Signal to Interference Ratio
%SINR = SNR/(INR + 1);          % Signal to Interference plus Noise Ratio


%% Define dB values and relative powers;
SNR_db = 10*log10(SNR);

P_x = 1;         % Normalize received Signal Power to 1
P_n = 1./SNR;    % noise power 
P_i = P_x/SIR;   % Interference Power (received)


%% Define Constellations
% QAM Constellations (normalized to average symbol power of 1)
N_Constellations = 2; % Considering 4QAM and 16QAM
Constellation_4  = [-1+1i; 1+1i; ...
                    -1-1i; 1-1i];
Constellation_4  = Constellation_4./sqrt(mean(abs(Constellation_4).^2));

Constellation_16 = [-3+3i; -1+3i; 1+3i; 3+3i; ...
                    -3+1i; -1+1i; 1+1i; 3+1i; ...
                    -3-1i; -1-1i; 1-1i; 3-1i; ...
                    -3-3i; -1-3i; 1-3i; 3-3i];
Constellation_16 = Constellation_16./sqrt(mean(abs(Constellation_16).^2));


%% Allocate Array
L_S_x   = zeros(Len_block,Num_blocks);
L_S_i   = zeros(Len_block,Num_blocks);
X       = zeros(Len_block,Num_blocks);
Y       = zeros(Len_block,Num_blocks);
I       = zeros(Len_block,Num_blocks);
I_rx    = zeros(Len_block,Num_blocks);
N       = zeros(Len_block,Num_blocks);

%% Generate signal and interference Symbols
% Block Settings (Contellation & presence of interference constant per block)
L_Constellations = randi(N_Constellations,1,Num_blocks);
L_Interference   = rand(1,Num_blocks) < p_int;

% Gen symbols for each block
for block = 1:Num_blocks

    % Define signal constellation for the block
    Constellation = L_Constellations(block);
    switch Constellation
        case 1
            Constellation_S = Constellation_4;
        case 2
            Constellation_S = Constellation_16;
        otherwise
            Constellation_S = Constellation_4;
    end
    
    % Define interference constellation for the block (QAM4 or none)
    Constellation_I = L_Interference(block)*Constellation_4;

    % Randomly select symbols for the given block
    L_S_x(:,block) = randi(length(Constellation_S),Len_block,1);
    L_S_i(:,block) = randi(length(Constellation_I),Len_block,1);

    X(:,block) = Constellation_S(L_S_x(:,block));
    I(:,block) = Constellation_I(L_S_i(:,block));

    %% Generate received signal
    I_rx(:,block) = sqrt(P_i)*I(:,block);                           % Received Interference
    Y(:,block)    = awgn(X(:,block),SNR_db(block)) + I_rx(:,block); % Noise dependent on SNR
    N(:,block)    = Y(:,block) - X(:,block) - I_rx(:,block);        % Noise values (Random)

    % Error Check received signal/interference/noise power
    P_n_actual(block) = mean(abs(N(:,block)).^2,'all');

end
% Error Check received signal/interference/noise power
P_i_actual = mean(abs(I_rx).^2,'all')/p_int;
P_x_actual = mean(abs(X).^2,'all');

if(PARAMS.NORMALIZE_RESULTS)
    % Normalize so that mean(abs(Y).^2,'all') = 1
    Y = Y./sqrt(mean(abs(Y).^2,'all'));
end




%% "Expert" Decoders
C_hat_4 = 1*ones(1,Num_blocks);
S_hat_4 = decode_QAM(Y,4);

C_hat_16 = 2*ones(1,Num_blocks);
S_hat_16 = decode_QAM(Y,16);

% Determine Error Rates for perfect modulation classification
SER_4  = 1 - sum(S_hat_4(:,C_hat_4==L_Constellations) == L_S_x(:,C_hat_4==L_Constellations),'all')./(Len_block*sum(C_hat_4==L_Constellations));
SER_16 = 1 - sum(S_hat_16(:,C_hat_16==L_Constellations) == L_S_x(:,C_hat_16==L_Constellations),'all')./(Len_block*sum(C_hat_16==L_Constellations));
       

%% Visualize
if PARAMS.VISUALIZE
    for i=1:min(5,Num_blocks)
        figure()
        scatter(real(Y(:,i)),imag(Y(:,i)))
        if(PARAMS.NORMALIZE_RESULTS)
            xlim([-1.5,1.5]);
            ylim([-1.5,1.5]);
        end
    end
    
    % Plot ideal constellations
    figure()
    scatter(real(Constellation_4),imag(Constellation_4),'d');
    hold on
    scatter(real(Constellation_16),imag(Constellation_16),'r');
    legend('4QAM','16QAM');

end

save my_labels.mat L_Constellations L_Interference L_S_x L_S_i
save my_decoder.mat S_hat_4 S_hat_16 SER_16 SER_4
save my_data.mat Y
