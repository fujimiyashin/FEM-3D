function res = unit_stiffness_F(X, Y, Z, force)
    res = zeros(24, 1);
    gauss_pnt = [-1/sqrt(3), 1/sqrt(3)];
    weight = [1, 1];
    for i=1:length(gauss_pnt)
        for j=1:length(gauss_pnt)
            for k=1:length(gauss_pnt)
                ksi = gauss_pnt(i);
                eta = gauss_pnt(j);
                zeta = gauss_pnt(k);
                jacobi = get_jacobi(ksi, eta, zeta, X, Y, Z);
                jacobi_det = det(jacobi);
                N1 = 0.125*(1-ksi)*(1-eta)*(1-zeta)*eye(3);
                N2 = 0.125*(1+ksi)*(1-eta)*(1-zeta)*eye(3);
                N3 = 0.125*(1+ksi)*(1+eta)*(1-zeta)*eye(3);
                N4 = 0.125*(1-ksi)*(1+eta)*(1-zeta)*eye(3);
                N5 = 0.125*(1-ksi)*(1-eta)*(1+zeta)*eye(3);
                N6 = 0.125*(1+ksi)*(1-eta)*(1+zeta)*eye(3);
                N7 = 0.125*(1+ksi)*(1+eta)*(1+zeta)*eye(3);
                N8 = 0.125*(1-ksi)*(1+eta)*(1+zeta)*eye(3);
                N = [N1, N2, N3, N4, N5, N6, N7, N8];

                res = res + N'*force*jacobi_det*weight(i)*weight(j)*weight(k);
            end
        end
    end
end

