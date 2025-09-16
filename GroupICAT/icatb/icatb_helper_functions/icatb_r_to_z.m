function z = icatb_r_to_z(r)
%function t = r_to_t(r)
%r=cor coeff

z = .5*[log(.00000001+1+r)-log(.00000001+1-r)];