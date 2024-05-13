% This function returns the equation for linear fits.
% 
% INPUT:
%   - bpar: coefficients (linear fit)
%   - xvar: variable (X-axis)
%   - yvar: variable (Y-axis)
%   - eqn_type: format ('text' or 'latex')
% 
% OUPUT:
%   - eqn_text: equation in selected format
% 
% Author: sebastien.viscardy@aeronomie.be
%%
function eqn_text= SS_mk_eqn_text(bpar,xvar,yvar,eqn_type)

%% Coefficients: y = ax + b
a = bpar(2);
b = bpar(1);

switch eqn_type
    case 'text'
        %% coefficient a
        if abs(log10(abs(a)))>=2
            if a >= 0
                atxt = SS_conv_expon4txt(a,2.2);
            else
                atxt = ['-',SS_conv_expon4txt(abs(a),2.2)];
            end
        else
            if a >= 0
                atxt = num2str(a,'%2.2f');
            else
                atxt = ['-',num2str(abs(a),'%2.2f')];
            end
        end
        
        %% coefficient b
        if abs(log10(abs(b)))>=3
            if b >= 0
                btxt = [' + ',SS_conv_expon4txt(b,2.2)];
            else
                btxt = [' - ',SS_conv_expon4txt(abs(b),2.2)];
            end
        else
            if b >= 0
                btxt = [' + ',num2str(b,'%4.2f')];
            else
                btxt = [' - ',num2str(abs(b),'%4.2f')];
            end
        end
        
        %% Make text of equation
        eqn_text = [yvar,' = ',atxt,' ',xvar,btxt];
        
    case 'latex'
        %% coefficient a
        if abs(log10(abs(a)))>=2
            if a >= 0
                atxt = SS_conv_expon4txt(a,2.2);
            else
                atxt = ['-',SS_conv_expon4txt(abs(a),2.2)];
            end
        else
            if a >= 0
                atxt = num2str(a,'%2.2f');
            else
                atxt = ['-',num2str(abs(a),'%2.2f')];
            end
        end
        
        %% coefficient b
        if abs(log10(abs(b)))>=3
            if b >= 0
                btxt = [' + ',SS_conv_expon4txt(b,2.2)];
            else
                btxt = [' - ',SS_conv_expon4txt(abs(b),2.2)];
            end
        else
            if b >= 0
                btxt = [' + ',num2str(b,'%4.2f')];
            else
                btxt = [' - ',num2str(abs(b),'%4.2f')];
            end
        end
        
        %% Make text of equation
        eqn_text = ['$$',yvar,' = ',atxt,' \ ',xvar,btxt,'$$'];
end