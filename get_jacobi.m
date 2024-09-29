function res = get_jacobi(ksi, eta, zeta, X, Y, Z)
    k = 1/8;
    jacobiA = k*[-(1-eta)*(1-zeta), (1-eta)*(1-zeta), (1+eta)*(1-zeta), -(1+eta)*(1-zeta), -(1-eta)*(1+zeta), (1-eta)*(1+zeta), (1+eta)*(1+zeta), -(1+eta)*(1+zeta);
                 -(1-ksi)*(1-zeta), -(1+ksi)*(1-zeta), (1+ksi)*(1-zeta), (1-ksi)*(1-zeta), -(1-ksi)*(1+zeta), -(1+ksi)*(1+zeta), (1+ksi)*(1+zeta), (1-ksi)*(1+zeta);
                 -(1-ksi)*(1-eta), -(1+ksi)*(1-eta), -(1+ksi)*(1+eta), -(1-ksi)*(1+eta), (1-ksi)*(1-eta), (1+ksi)*(1-eta), (1+ksi)*(1+eta), (1-ksi)*(1+eta);];
    jacobiB = [X, Y, Z];
    res = jacobiA*jacobiB;
end

