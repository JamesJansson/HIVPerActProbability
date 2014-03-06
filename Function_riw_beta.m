function l = Function_riw_beta(v, Params)
global n y n_ins_neg n_ins_unknown n_ins_pos n_recw_neg n_recw_unknown n_recw_pos n_rec_neg n_rec_unknown n_rec_pos c lastones undiag_u undiag_n;


if Params.TrackingCircumcision==1
    beta_r=v(1);
    beta_i=v(2);
    beta_w=v(3);
    beta_ic=v(4);
elseif Params.TrackingWithdrawal==1
    beta_r=v(1);
    beta_i=v(2);
    beta_w=v(3);
    beta_ic=v(2);
else
    beta_r=v(1);
    beta_i=v(2);
    beta_w=v(1);
    beta_ic=v(2);
end


temp_i = 0;
temp_u = 0;
f=zeros(n,1);


i=y==1;%infected
fi=log(1-(1-beta_r).^(n_rec_pos(i,1)+ undiag_u*n_rec_unknown(i,1)+ undiag_n*n_rec_neg(i,1)).*(1-beta_w).^(n_recw_pos(i,1)+ undiag_u*n_recw_unknown(i,1)+ undiag_n*n_recw_neg(i,1)).*(1-(1-c(i,1))*beta_i-beta_ic*c(i,1)).^(n_ins_pos(i,1)+undiag_u*n_ins_unknown(i,1)+undiag_n*n_ins_neg(i,1)));
i=y==0;%uninfected
fu=log(1-beta_r).*(n_rec_pos(i,1)+undiag_u*n_rec_unknown(i,1)+undiag_n*n_rec_neg(i,1))+log(1-beta_w).*(n_recw_pos(i,1)+undiag_u*n_recw_unknown(i,1)+undiag_n*n_recw_neg(i,1))+log(1-(1-c(i,1)).*beta_i-beta_ic*c(i,1)).*(n_ins_pos(i,1)+undiag_u*n_ins_unknown(i,1)+undiag_n*n_ins_neg(i,1));

l=-(sum(fi)+sum(fu));


% for i=1:n
% 
%     %DUAL (receptive-insertive) BETA CASE
%     
%     if y(i,1)==1 %if HIV positive 
%         
%         f(i)=1-(1-beta_r).^(n_rec_pos(i,1)+ undiag_u*n_rec_unknown(i,1)+ undiag_n*n_rec_neg(i,1))*(1-beta_w).^(n_recw_pos(i,1)+ undiag_u*n_recw_unknown(i,1)+ undiag_n*n_recw_neg(i,1))*(1-(1-c(i,1))*beta_i-beta_ic*c(i,1)).^(n_ins_pos(i,1)+undiag_u*n_ins_unknown(i,1)+undiag_n*n_ins_neg(i,1));
%         
%         if f(i)>0 %to avoid log 0 errors where people get HIV without sexual contact with a known positive partner
%         
%             temp_i=temp_i+log(f(i));
% 
%         end
%     else% if HIV negaitive
% 
%         
%         temp_u=temp_u+log(1-beta_r)*(n_rec_pos(i,1)+undiag_u*n_rec_unknown(i,1)+undiag_n*n_rec_neg(i,1))+log(1-beta_w)*(n_recw_pos(i,1)+undiag_u*n_recw_unknown(i,1)+undiag_n*n_recw_neg(i,1))+log(1-(1-c(i,1))*beta_i-beta_ic*c(i,1))*(n_ins_pos(i,1)+undiag_u*n_ins_unknown(i,1)+undiag_n*n_ins_neg(i,1));
%     end
%     
% end
% 
% 
% l=-1*(temp_u+temp_i); %-1 times because we want to maximize 'l' and matlab minimizes functions 
