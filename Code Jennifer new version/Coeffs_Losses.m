function [C_rp , C_rq, C_xp, C_xq] =Coeffs_Losses(YY, E_LF, COEFFcomplex_alph, COEFFqcomplex_alph, nodesnoslack)

% nodesnoslack has to be a vector of all the nodes of the system.. in reality if we know 
% that we computed all the SC in all nodes we can omit this variable and just write
% l = 1:length(COEFFcomplex_alph)

temp_p = zeros(1,length(nodesnoslack));
temp_q = zeros(1,length(nodesnoslack));

% test = 0;
% for i = 1:size(E_LF,1)
%     tmp_YijEj = 0;
%     for j = 1:size(E_LF,1)
%         tmp_YijEj = tmp_YijEj + conj(YY(i,j))*conj(E_LF(j));
%     end
%     test = test + E_LF(i)*tmp_YijEj;
% end

for l = 1:length(nodesnoslack)
    
    for i = 1:size(E_LF,1)
        
        tmp_YijEj = 0;
        tmp_YijdEjdP = 0;
        tmp_YijdEjdQ = 0;
        
        for j = 1:size(E_LF,1)
            tmp_YijEj = tmp_YijEj + conj(YY(i,j))*conj(E_LF(j));
            tmp_YijdEjdP = tmp_YijdEjdP + conj(YY(i,j))*conj(COEFFcomplex_alph(j,l));
            tmp_YijdEjdQ = tmp_YijdEjdQ + conj(YY(i,j))*conj(COEFFqcomplex_alph(j,l));
        end 
        
        temp_p(l) = temp_p(l) + COEFFcomplex_alph(i,l)*tmp_YijEj + E_LF(i)*tmp_YijdEjdP;   
        temp_q(l) = temp_q(l) + COEFFqcomplex_alph(i,l)*tmp_YijEj + E_LF(i)*tmp_YijdEjdQ;
    end 
end

C_rp = -real(temp_p); 
C_rq = -real(temp_q);
C_xp = -imag(temp_p);
C_xq = -imag(temp_q);

end