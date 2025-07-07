 function [results,R2]=cosinefit2(rates,directions,weights)
%rates is a 8 numbers vector of firing rate for each of the 8
%direction(order same as in "directions below-0,45,90...)

%this fiunction 1. create a cosine model
model=fittype('a*cos(x*pi/180+b)+c');
%directions=[0 45 90 135 180 225 270 315];
opts = fitoptions;
opts.weights=weights;
%2. find R of correlation between data and model
[results, gof]=fit(directions',rates',model,opts);
% R2=gof.adjrsquare;
R2=gof.rsquare;
% if R2<0
%     'wait'
% end
%notes: R is found by methods of non linear least squares. to change use
%fitoptions function
%second R@ here is the adjusted R2, i.e it can take negative values. if you
%the non-adjusted-always -positive method, use R@=gof.rsquare
