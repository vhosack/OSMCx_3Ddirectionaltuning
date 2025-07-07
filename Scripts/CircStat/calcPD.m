function [BD, AMP,BASE, R2]=calcPD(rates,directions,method)
%function [BD, AMP,BASE, R2]=computePD(rates,directions,method)
% directions should be in angles when 0 is the X axis  to the right
% rates is a cell array with length as directions and each cell is an array of
% trial rates
% or rates can be an array of the mean rates in each direction

MIN_NUM_OF_SUCC_TRIALS=3;
MIN_NUM_OF_DIRECTIONS=4;

if length(rates) ~= length(directions)
    warning('directions and rates have different length');
    PD_theta=0;
    PD_size=0;
    return;
end

if iscell(rates) % cell array of trials rates
    % remove directions with too few trials
    d=1;
    for i=1:length(rates)
        x=find(~isnan(rates{i}));
        if length(x) >= MIN_NUM_OF_SUCC_TRIALS
            RATES{d}=rates{i}(x);
            DIRECTIONS(d)=directions(i);
            d=d+1;
        end
    end
    for i=1:length(RATES)
        M_RATES(i) = mean(RATES{i});
    end
else
    x=find(~isnan(rates));
    M_RATES=rates(x);
    DIRECTIONS=directions(x);
end

% abort if not enough directions
if length(~isnan(M_RATES)) < MIN_NUM_OF_DIRECTIONS
    BD = -1; AMP=0; BASE=0; R2=0;
    warning('too few directions for calculating PD');
    return;
end
       
% calculate PD according to COS fit or vector summation
if nargin > 2 & strcmp(lower(method),'cos')
%     [BD, AMP,BASE, R2] = fitcos(M_RATES, DIRECTIONS,'plotres');
    [BD, AMP,BASE, R2] = fitcos(M_RATES, DIRECTIONS,[]);
%     [results,R2] = cosinefit2(M_RATES, DIRECTIONS,[]);
    
else
    R2=1;

    dir_cos=cos(2*pi*DIRECTIONS/360);
    dir_sin=sin(2*pi*DIRECTIONS/360);
    Xs=M_RATES.*dir_cos;
    Ys=M_RATES.*dir_sin;

    AMP=sqrt(sum(Ys)^2+sum(Xs)^2);
    BD=atan2(sum(Ys),sum(Xs));
    if BD<0
        BD= 2*pi + BD;
    end
    AMP=AMP/length(DIRECTIONS);
    BASE=mean(M_RATES);
end

BD=360*BD/2/pi;

if BD<0 BD=360+BD; end

        

