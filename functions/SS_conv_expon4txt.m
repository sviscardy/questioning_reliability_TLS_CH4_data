% val_txt = SS_conv_expon4txt(val0,format)
%
% This function returns the version for text of an exponential value in the
% format chosen.
%
% Ex.: val0 = 1.432e-15
% val_txt = SS_conv_expon4txt(val0,2.2) ==> '1.43\times10^{-15}'
%
% INPUT:
%   - val0: initial exponential value
%   - format: format of argument
%       (ex.: 2.2 --> 1.432e-15 ==> 1.43\times10^{-15})
%
% OUTPUT:
%   - exponential value (string) that can be used for text.
%
% Author: sebastien.viscardy@aeronomie.be
%%
function val_txt = SS_conv_expon4txt(val0,format)

sign_val = abs(val0)/val0;
val0     = abs(val0);
exp_val  = floor(log10(val0));
arg_val  = val0/10^exp_val;

if ( sign_val == -1 )
    val_txt = ['-',num2str(arg_val,['%',num2str(format),'f']), ...
        '\times10^{',num2str(exp_val),'}'];
elseif ( sign_val == 1 )
    val_txt = [num2str(arg_val,['%',num2str(format),'f']), ...
        '\times10^{',num2str(exp_val),'}'];
end
end