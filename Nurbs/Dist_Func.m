function [gg,dgg]=Dist_Func(h,gpt,imix)

switch imix
    case 1
        gg = gpt-4/3*gpt^3/h^2;         %%% Model Reddy
        dgg= 1-4*gpt^2/h^2;             %%% Model Reddy
    case 2
        gg = sin(pi*gpt/h);             %%% Model Arya
        dgg= pi/h*cos(pi*gpt/h);        %%% Model Arya
    case 3
        gg = (h/pi)*sin(pi*gpt/h);      %%% Model Touratier
        dgg= cos(pi*gpt/h);             %%% Model Touratier
    case 4
        gg = gpt*exp(-2*(gpt/h)^2);             %%% Model Karama
        dgg= (1-4*(gpt/h)^2)*exp(-2*(gpt/h)^2); %%% Model Karama
    case 5
        gg = -h*sinh(gpt/h)+gpt*cosh(1/2);   %%% Model Soldatos
        dgg= -cosh(gpt/h)+cosh(1/2);         %%% Model Soldatos
    case 6
        alpha=3;
        somu =-2*(gpt/h)^2/log(alpha);
        gg   = gpt*alpha^somu;               %%% Model Aydogdu
        dgg  = alpha^somu*(1-4*(gpt/h)^2);   %%% Model Aydogdu
    case 7
        gg = h*atan(2*gpt/h)-gpt;               % % chien  atan
        dgg= 2/(1+(2*gpt/h)^2)-1;
    case 8
        gg = atan(sin(pi*gpt/h));               % % chien mixed atan
        dgg= pi/h*cos(pi*gpt/h)/(1+(sin(pi*gpt/h))^2);
    case 9
        gg = asinh(sin(pi*gpt/h));              % % chien sin hyperbolic
        dgg= pi/h*cos(pi*gpt/h)/sqrt(1+(sin(pi*gpt/h))^2);
    case 10
        gg = 7/8*gpt-2/h^2*gpt^3+2/h^4*gpt^5;             % % hung bac 5
        dgg= 7/8-6/h^2*gpt^2+10/h^4*gpt^4;
end