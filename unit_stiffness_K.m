function res = unit_stiffness_K(G, nu, X, Y, Z)
res = zeros(24, 24);
E = G*2*(1+nu);
lame = E*nu/((1+nu)*(1-2*nu));
D1 = [lame+2*G, lame, lame;
      lame, lame+2*G, lame;
      lame, lame, lame+2*G];
D2 = zeros(3,3);
D3 = G*eye(3);
D = [D1, D2;
    D2, D3];
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
            N_ksi_eta_zeta = 0.125*[-(1-eta)*(1-zeta), (1-eta)*(1-zeta), (1+eta)*(1-zeta), -(1+eta)*(1-zeta), -(1-eta)*(1+zeta), (1-eta)*(1+zeta), (1+eta)*(1+zeta), -(1+eta)*(1+zeta);% 对ksi求导
                                    -(1-ksi)*(1-zeta), -(1+ksi)*(1-zeta), (1+ksi)*(1-zeta), (1-ksi)*(1-zeta), -(1-ksi)*(1+zeta), -(1+ksi)*(1+zeta), (1+ksi)*(1+zeta), (1-ksi)*(1+zeta);% 对eta求导
                                    -(1-ksi)*(1-eta), -(1+ksi)*(1-eta), -(1+ksi)*(1+eta), -(1-ksi)*(1+eta), (1-ksi)*(1-eta), (1+ksi)*(1-eta), (1+ksi)*(1+eta), (1-ksi)*(1+eta);];% 对zeta求导
            N_ksi_eta_zeta = N_ksi_eta_zeta';
            N1_xyz = jacobi\[N_ksi_eta_zeta(1,1);N_ksi_eta_zeta(1,2);N_ksi_eta_zeta(1,3)];
            N2_xyz = jacobi\[N_ksi_eta_zeta(2,1);N_ksi_eta_zeta(2,2);N_ksi_eta_zeta(2,3)];
            N3_xyz = jacobi\[N_ksi_eta_zeta(3,1);N_ksi_eta_zeta(3,2);N_ksi_eta_zeta(3,3)];
            N4_xyz = jacobi\[N_ksi_eta_zeta(4,1);N_ksi_eta_zeta(4,2);N_ksi_eta_zeta(4,3)];
            N5_xyz = jacobi\[N_ksi_eta_zeta(5,1);N_ksi_eta_zeta(5,2);N_ksi_eta_zeta(5,3)];
            N6_xyz = jacobi\[N_ksi_eta_zeta(6,1);N_ksi_eta_zeta(6,2);N_ksi_eta_zeta(6,3)];
            N7_xyz = jacobi\[N_ksi_eta_zeta(7,1);N_ksi_eta_zeta(7,2);N_ksi_eta_zeta(7,3)];
            N8_xyz = jacobi\[N_ksi_eta_zeta(8,1);N_ksi_eta_zeta(8,2);N_ksi_eta_zeta(8,3)];

            B1 = [N1_xyz(1), 0, 0;
                  0, N1_xyz(2), 0;
                  0, 0, N1_xyz(3);
                  N1_xyz(2), N1_xyz(1), 0;
                  0, N1_xyz(3), N1_xyz(2);
                  N1_xyz(3), 0, N1_xyz(1)];

            B2 = [N2_xyz(1), 0, 0;
                  0, N2_xyz(2), 0;
                  0, 0, N2_xyz(3);
                  N2_xyz(2), N2_xyz(1), 0;
                  0, N2_xyz(3), N2_xyz(2);
                  N2_xyz(3), 0, N2_xyz(1)];

            B3 = [N3_xyz(1), 0, 0;
                  0, N3_xyz(2), 0;
                  0, 0, N3_xyz(3);
                  N3_xyz(2), N3_xyz(1), 0;
                  0, N3_xyz(3), N3_xyz(2);
                  N3_xyz(3), 0, N3_xyz(1)];

            B4 = [N4_xyz(1), 0, 0;
                  0, N4_xyz(2), 0;
                  0, 0, N4_xyz(3);
                  N4_xyz(2), N4_xyz(1), 0;
                  0, N4_xyz(3), N4_xyz(2);
                  N4_xyz(3), 0, N4_xyz(1)];

            B5 = [N5_xyz(1), 0, 0;
                  0, N5_xyz(2), 0;
                  0, 0, N5_xyz(3);
                  N5_xyz(2), N5_xyz(1), 0;
                  0, N5_xyz(3), N5_xyz(2);
                  N5_xyz(3), 0, N5_xyz(1)];

            B6 = [N6_xyz(1), 0, 0;
                  0, N6_xyz(2), 0;
                  0, 0, N6_xyz(3);
                  N6_xyz(2), N6_xyz(1), 0;
                  0, N6_xyz(3), N6_xyz(2);
                  N6_xyz(3), 0, N6_xyz(1)];

            B7 = [N7_xyz(1), 0, 0;
                  0, N7_xyz(2), 0;
                  0, 0, N7_xyz(3);
                  N7_xyz(2), N7_xyz(1), 0;
                  0, N7_xyz(3), N7_xyz(2);
                  N7_xyz(3), 0, N7_xyz(1)];

            B8 = [N8_xyz(1), 0, 0;
                  0, N8_xyz(2), 0;
                  0, 0, N8_xyz(3);
                  N8_xyz(2), N8_xyz(1), 0;
                  0, N8_xyz(3), N8_xyz(2);
                  N8_xyz(3), 0, N8_xyz(1)];

            B = [B1, B2, B3, B4, B5, B6, B7, B8];

            res = res + B'*D*B*jacobi_det*weight(i)*weight(j)*weight(k);
        end
    end
end
end

