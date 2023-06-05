% REDO EASIER 3 for loops!

function [Icoeff, Icoeff_complex, Icoeffq,Icoeffq_complex, Currents]=Coeffs_Currents(YYL,YYT,E,COEFFcomplex_alph,COEFFqcomplex_alph,nph)
% !!!! YYT = A'*Yl*A
    Currents = cell(size(YYL)/nph);
    Icoeff = cell(1,size(COEFFcomplex_alph,2));
    Icoeffq = cell(1,size(COEFFcomplex_alph,2));
    Icoeff_complex = cell(1,size(COEFFcomplex_alph,2));
    Icoeffq_complex = cell(1,size(COEFFcomplex_alph,2));
    
    % For all active/reactive power injections
    for l = 1:size(COEFFcomplex_alph,2)
        % Re-initialize Temp Matrices
        temp_p = cell(size(YYL)/nph);
        temp_q = cell(size(YYL)/nph); 
        temp_p_complex = cell(size(YYL)/nph);
        temp_q_complex = cell(size(YYL)/nph);
        
        % For all rows of Admittance Matrix
        for i = 1:nph:size(YYL,1)
            for j = 1:nph:size(YYL,2)
                
                % Compute Currents only once
                if (l == 1)
                    Currents{ceil(i/nph),ceil(j/nph)} = YYL(i:i+(nph-1),j:j+(nph-1))*(E(i:i+(nph-1)) - E(j:j+(nph-1))) + ...
                                    YYT(i:i+(nph-1),j:j+(nph-1))*E(i:i+(nph-1));
                end
                
                if( max(abs(Currents{ceil(i/nph),ceil(j/nph)})) ~= 0)
                    % Compute Complex Coeffs
                    dIij_dPl = YYL(i:i+(nph-1),j:j+(nph-1))*(COEFFcomplex_alph(i:i+(nph-1),l) - COEFFcomplex_alph(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))*COEFFcomplex_alph(i:i+(nph-1),l);
                    dIij_dQl = YYL(i:i+(nph-1),j:j+(nph-1))*(COEFFqcomplex_alph(i:i+(nph-1),l) - COEFFqcomplex_alph(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))*COEFFqcomplex_alph(i:i+(nph-1),l);

                    % Compute Coeffs
                    temp_p{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dPl));
                    temp_q{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dQl));
                    temp_p_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dPl;
                    temp_q_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dQl;
                else
                    temp_p{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_q{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_p_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_q_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                end
                
            end
        end
        
        Icoeff{l} = temp_p;
        Icoeffq{l} = temp_q;
        Icoeff_complex{l} = temp_p_complex;
        Icoeffq_complex{l} = temp_q_complex; 
    end
        
end