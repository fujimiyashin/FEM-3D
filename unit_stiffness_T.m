function res = unit_stiffness_T(ksi_lis, eta_lis, zeta_lis, weight, X, Y, Z, force)
    res = zeros(24, 1);
    for i=1:length(ksi_lis)
        for j=1:length(eta_lis)
            for k=1:length(zeta_lis)
                ksi = ksi_lis(i);
                eta = eta_lis(j);
                zeta = zeta_lis(k);

                N1 = 0.125*(1-ksi)*(1-eta)*(1-zeta)*eye(3);
                N2 = 0.125*(1+ksi)*(1-eta)*(1-zeta)*eye(3);
                N3 = 0.125*(1+ksi)*(1+eta)*(1-zeta)*eye(3);
                N4 = 0.125*(1-ksi)*(1+eta)*(1-zeta)*eye(3);
                N5 = 0.125*(1-ksi)*(1-eta)*(1+zeta)*eye(3);
                N6 = 0.125*(1+ksi)*(1-eta)*(1+zeta)*eye(3);
                N7 = 0.125*(1+ksi)*(1+eta)*(1+zeta)*eye(3);
                N8 = 0.125*(1-ksi)*(1+eta)*(1+zeta)*eye(3);
                N = [N1, N2, N3, N4, N5, N6, N7, N8];
                
                if abs(ksi)==1
                    jacobi = get_jacobi(ksi, eta, zeta, X, Y, Z);
                    a = jacobi(2,2)*jacobi(3,3)-jacobi(2,3)*jacobi(3,2);
                    b = jacobi(2,1)*jacobi(3,3)-jacobi(2,3)*jacobi(3,1);
                    c = jacobi(2,1)*jacobi(3,2)-jacobi(2,2)*jacobi(3,1);
                    jacobi_det = sqrt(power(a, 2)+power(b, 2)+power(c, 2));
                    res = res + N'*force*jacobi_det*weight(j)*weight(k);
                elseif abs(eta)==1
                    jacobi = get_jacobi(ksi, eta, zeta, X, Y, Z);
                    a = jacobi(1,2)*jacobi(3,3)-jacobi(1,3)*jacobi(3,2);
                    b = jacobi(1,1)*jacobi(3,3)-jacobi(1,3)*jacobi(3,1);
                    c = jacobi(1,1)*jacobi(3,2)-jacobi(1,2)*jacobi(3,1);
                    jacobi_det = sqrt(power(a, 2)+power(b, 2)+power(c, 2));
                    res = res + N'*force*jacobi_det*weight(i)*weight(k);
                elseif abs(zeta)==1
                    jacobi = get_jacobi(ksi, eta, zeta, X, Y, Z);
                    a = jacobi(1,2)*jacobi(2,3)-jacobi(1,3)*jacobi(2,2);
                    b = jacobi(1,1)*jacobi(2,3)-jacobi(1,3)*jacobi(2,1);
                    c = jacobi(1,1)*jacobi(2,2)-jacobi(1,2)*jacobi(2,1);
                    jacobi_det = sqrt(power(a, 2)+power(b, 2)+power(c, 2));
                    res = res + N'*force*jacobi_det*weight(j)*weight(i);
                end
            end
        end
    end
end

