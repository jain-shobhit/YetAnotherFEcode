module DpROM

using Dates
using TensorOperations
using LinearAlgebra
using SparseArrays

export red_stiff_tensors, red_stiff_tensors_defects
export stiffness_matrix_derivative, stiffness_matrix_sensitivity
export stiffness_matrix_linear

# ASSEMBLY FUNCTIONS ___________________________________________________________
# ASSEMBLE ROM stiffness tensors (standard, no parametrization)
function red_stiff_tensors(elements, nodes, connectivity, C, V, XGauss, WGauss)
    time1 = now()

    # useful dimensions
    nD = size(nodes,2)
    nel = size(elements,1)
    ndofs = size(connectivity,2)
    nw = length(WGauss)
    nv = size(V,2)

    G_FUN = Gselector(nD,ndofs)

    if nD==3
        L1, L2, L31, H = constantMatrices3D()
    elseif nD==2
        L1, L2, L31, H = constantMatrices2D()
    end

    # initialize tensors
    K2n  = zeros(nv,nv)
    K3n  = zeros(nv,nv,nv)
    K4n  = zeros(nv,nv,nv,nv)

    # cycle over all the elements
    for e in 1:nel
        el_nodes = elements[e,:]        # IDs of the element nodes
        el_dofs  = connectivity[e,:]    # IDs of the element dofs
        xyz = nodes[el_nodes,:]         # element's coordinates x, y, z
        Ve = V[el_dofs,:]               # Projection Basis V (only element's dofs)

        # cycle over Gauss quadrature points (integration over volume)
        for i in 1:nw
            ISO = XGauss[:,i]           # Gauss point coordinates
            w = WGauss[i]               # Gauss weight
            G, detJ = G_FUN(ISO,xyz)    # shape function derivative matrix G and Jacobian
            # compute the element-level reduced tensors at the Gauss point
            CwJ = C .* (w * detJ);
            K2nE, K3nE, K4nE = element_tensor_standard(G,CwJ,H,L1,Ve)
            # sum over the element contributions
            K2n .+= K2nE
            K3n .+= K3nE
            K4n .+= K4nE
        end
    end
    time2 = now()
    time3 = time2 - time1
    totaltime = time3.value
    # return all the tensors
    K2n, K3n, K4n, totaltime;
end
# ASSEMBLE DpROM stiffness tensors
function red_stiff_tensors_defects(formulation, volume, elements, nodes, connectivity, C, V, U, XGauss, WGauss)
    time1 = now()

    # useful dimensions
    nD = size(nodes,2)
    nel = size(elements,1)
    ndofs = size(connectivity,2)
    nw = length(WGauss)
    nv = size(V,2)
    nu = size(U,2)

    G_FUN = Gselector(nD,ndofs)

    element_tensor = element_tensor_selector(formulation)

    if nD==3
        L1, L2, L31, H = constantMatrices3D()
        compute_defect_divergence = compute_defect_divergence_3D
    elseif nD==2
        L1, L2, L31, H = constantMatrices2D()
        compute_defect_divergence = compute_defect_divergence_2D
    end

    # initialize tensors
    K2n  = [zeros(nv,nv)             for _ in 1:1, _ in 1:(nu+1)]
    K3d  = [zeros(nv,nv,nu)          for _ in 1:1, _ in 1:(nu+1)]
    K4dd = [zeros(nv,nv,nu,nu)       for _ in 1:1, _ in 1:(nu+1)]
    K3n  = [zeros(nv,nv,nv)          for _ in 1:1, _ in 1:(nu+1)]
    K4d  = [zeros(nv,nv,nu,nv)       for _ in 1:1, _ in 1:(nu+1)]
    K5dd = [zeros(nv,nv,nu,nv,nu)    for _ in 1:1, _ in 1:(nu+1)]
    K4n  = [zeros(nv,nv,nv,nv)       for _ in 1:1, _ in 1:(nu+1)]
    K5d  = [zeros(nv,nv,nv,nu,nv)    for _ in 1:1, _ in 1:(nu+1)]
    K6dd = [zeros(nv,nv,nu,nv,nu,nv) for _ in 1:1, _ in 1:(nu+1)]

    # cycle over all the elements
    for e in 1:nel
        el_nodes = elements[e,:]        # IDs of the element nodes
        el_dofs  = connectivity[e,:]    # IDs of the element dofs
        xyz = nodes[el_nodes,:]         # element's coordinates x, y, z
        Ve = V[el_dofs,:]               # Projection Basis V (only element's dofs)
        Ue = U[el_dofs,:]               # Defect Basis V (only element's dofs)

        # cycle over Gauss quadrature points (integration over volume)
        for i in 1:nw
            ISO = XGauss[:,i]           # Gauss point coordinates
            w = WGauss[i]               # Gauss weight
            G, detJ = G_FUN(ISO,xyz)    # shape function derivative matrix G and Jacobian
            CwJ = C .* (w * detJ);
            # compute the element-level reduced tensors at the Gauss point
            K2nE, K3dE, K4ddE, K3nE, K4dE, K5ddE, K4nE, K5dE, K6ddE, = element_tensor(G,CwJ,H,L1,L2,L31,Ve,Ue)
            # sum over the element contributions
            K2n[1]  .+= K2nE
            K3d[1]  .+= K3dE
            K4dd[1] .+= K4ddE
            K3n[1]  .+= K3nE
            K4d[1]  .+= K4dE
            K5dd[1] .+= K5ddE
            K4n[1]  .+= K4nE
            K5d[1]  .+= K5dE
            K6dd[1] .+= K6ddE

            # compute additional contributions (1 matrix for each defect)
            if volume==1
                for dd in 2 : (nu+1)
                    def_div = compute_defect_divergence(G, Ue[:,dd-1]);
                    K2n[dd]  .+= K2nE  .* def_div
                    K3d[dd]  .+= K3dE  .* def_div
                    K4dd[dd] .+= K4ddE .* def_div
                    K3n[dd]  .+= K3nE  .* def_div
                    K4d[dd]  .+= K4dE  .* def_div
                    K5dd[dd] .+= K5ddE .* def_div
                    K4n[dd]  .+= K4nE  .* def_div
                    K5d[dd]  .+= K5dE  .* def_div
                    K6dd[dd] .+= K6ddE .* def_div
                end
            end
        end
    end
    time2 = now()
    time3 = time2 - time1
    totaltime = time3.value
    # return all the tensors
    K2n, K3d, K4dd, K3n, K4d, K5dd, K4n, K5d, K6dd, totaltime;
end

# approiximated determinant (for integration over the defected volume)
function compute_defect_divergence_2D(G,Ue)
    TH = G*Ue;
    detF1a = TH[1] + TH[4]
    detF1a
end
function compute_defect_divergence_3D(G,Ue)
    TH = G*Ue;
    detF1a = TH[1] + TH[5] + TH[9]
    detF1a
end


# TENSOR FORMULATIONS __________________________________________________________
function element_tensor_selector(formulation)
    if formulation == "N1"
        element_tensor = element_tensor_N1
    elseif formulation == "N1T"
        element_tensor = element_tensor_N1t
    elseif formulation == "N0"
        element_tensor = element_tensor_N0
    elseif formulation == "STANDARD"
        element_tensor = element_tensor_standard
    end
    element_tensor
end
# NEUMANN EXPANSION order 1
# Reduced Tensors at element-level (evaluated at one gauss point)
function element_tensor_N0(G,C,H,L1,L2,L31,V,U)
    GHC = G'*H'*C
    CHG = C*H*G;
    @tensoropt begin
        # second order tensor K2
        K2n[I,J] := V[II,I]*G[i,II]*H[j,i]*C[j,k]*H[k,l]*G[l,JJ]*V[JJ,J];
        K3d[I,J,K] := V[II,I]*G[i,II]*(H[j,i]*C[j,k]*L1[k,l,a]*G[a,KK]*U[KK,K] + L1[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*H[k,l])*G[l,JJ]*V[JJ,J];
        K4dd[I,J,K,L] := V[II,I]*G[i,II]*L1[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*L1[k,l,b]*G[b,LL]*U[LL,L]*G[l,JJ]*V[JJ,J];

        # third order tensor K3
        K3n[I,J,K] := 0.5*V[II,I]*GHC[II,i]*L1[i,j,l]*G[l,KK]*V[KK,K]*G[j,JJ]*V[JJ,J] + V[II,I]*G[i,II]*L1[j,i,k]*G[k,KK]*V[KK,K]*CHG[j,JJ]*V[JJ,J];
        K4d[I,J,K,L] := V[II,I]*G[i,II]*( 0.5*L1[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*L1[k,l,b]*G[b,LL]*V[LL,L] + L1[j,i,a]*G[a,LL]*V[LL,L]*C[j,k]*L1[k,l,b]*G[b,KK]*U[KK,K] )*G[l,JJ]*V[JJ,J];

        # fourth order tensor K4
        K4n[I,J,K,L] := 0.5*V[II,I]*G[i,II]*L1[j,i,k]*G[k,JJ]*V[JJ,J] *C[j,l]* L1[l,m,n]*G[n,KK]*V[KK,K]*G[m,LL]*V[LL,L];
    end
    nv = size(V,2)
    nu = size(U,2)
    K5dd = zeros(nv,nv,nu,nv,nu)
    K5d  = zeros(nv,nv,nv,nu,nv)
    K6dd = zeros(nv,nv,nu,nv,nu,nv)
    K2n, K3d, K4dd, K3n, K4d, K5dd, K4n, K5d, K6dd
end
# NEUMAN EXPANSION order 1, truncated version
# Reduced Tensors at element-level (evaluated at one gauss point)
function element_tensor_N1t(G,C,H,L1,L2,L31,V,U)
    GHC = G'*H'*C
    CHG = C*H*G;
    @tensoropt begin
        # second order tensor K2
        K2n[I,J] := V[II,I]*G[i,II]*H[j,i]*C[j,k]*H[k,l]*G[l,JJ]*V[JJ,J];
        K3d[I,J,K] := V[II,I]*G[i,II]*(H[j,i]*C[j,k]*L2[k,l,a]*G[a,KK]*U[KK,K] + L2[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*H[k,l])*G[l,JJ]*V[JJ,J];
        K4dd[I,J,K,L] := V[II,I]*G[i,II]*L2[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*L2[k,l,b]*G[b,LL]*U[LL,L]*G[l,JJ]*V[JJ,J];

        # third order tensor K3
        K3n[I,J,K] := 0.5*V[II,I]*GHC[II,i]*L1[i,j,l]*G[l,KK]*V[KK,K]*G[j,JJ]*V[JJ,J] + V[II,I]*G[i,II]*L1[j,i,k]*G[k,KK]*V[KK,K]*CHG[j,JJ]*V[JJ,J];
        K4d[I,J,K,L] := V[II,I]*G[i,II]*( 0.5*L2[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*L1[k,l,b]*G[b,LL]*V[LL,L] + L1[j,i,a]*G[a,LL]*V[LL,L]*C[j,k]*L2[k,l,b]*G[b,KK]*U[KK,K] )*G[l,JJ]*V[JJ,J];

        # fourth order tensor K4
        K4n[I,J,K,L] := 0.5*V[II,I]*G[i,II]*L1[j,i,k]*G[k,JJ]*V[JJ,J] *C[j,l]* L1[l,m,n]*G[n,KK]*V[KK,K]*G[m,LL]*V[LL,L];
    end
    nv = size(V,2)
    nu = size(U,2)
    K5dd = zeros(nv,nv,nu,nv,nu)
    K5d  = zeros(nv,nv,nv,nu,nv)
    K6dd = zeros(nv,nv,nu,nv,nu,nv)
    K2n, K3d, K4dd, K3n, K4d, K5dd, K4n, K5d, K6dd
end
# NEUMANN EXPANSION order 0 (Budiansky, if integrated on nominal volume)
# Reduced Tensors at element-level (3d, evaluated at one gauss point)
function element_tensor_N1(G,C,H,L1,L2,L31,V,U)
    GHC = G'*H'*C
    CHG = C*H*G;
    @tensoropt begin
        # second order tensor K2
        K2n[I,J] := V[II,I]*G[i,II]*H[j,i]*C[j,k]*H[k,l]*G[l,JJ]*V[JJ,J];
        K3d[I,J,K] := V[II,I]*G[i,II]*(H[j,i]*C[j,k]*L2[k,l,a]*G[a,KK]*U[KK,K] + L2[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*H[k,l])*G[l,JJ]*V[JJ,J];
        K4dd[I,J,K,L] := V[II,I]*G[i,II]*L2[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*L2[k,l,b]*G[b,LL]*U[LL,L]*G[l,JJ]*V[JJ,J];

        # third order tensor K3
        K3n[I,J,K] := 0.5*V[II,I]*GHC[II,i]*L1[i,j,l]*G[l,KK]*V[KK,K]*G[j,JJ]*V[JJ,J] + V[II,I]*G[i,II]*L1[j,i,k]*G[k,KK]*V[KK,K]*CHG[j,JJ]*V[JJ,J];
        K4d[I,J,K,L] := V[II,I]*G[i,II]*(0.5*L2[j,i,a]*G[a,KK]*U[KK,K]*C[j,k]*L1[k,l,b]*G[b,LL]*V[LL,L] + L1[j,i,a]*G[a,LL]*V[LL,L]*C[j,k]*L2[k,l,b]*G[b,KK]*U[KK,K] + 2.0*L31[j,i,a,b]*G[b,KK]*U[KK,K]*G[a,LL]*V[LL,L]*C[j,k]*H[k,l] + H[j,i]*C[j,k]*L31[k,l,a,b]*G[b,KK]*U[KK,K]*G[a,LL]*V[LL,L])*G[l,JJ]*V[JJ,J];
        K5dd[I,J,K,L,M] := V[II,I]*G[i,II]*(2.0*L31[j,i,a,b]*G[b,KK]*U[KK,K]*G[a,LL]*V[LL,L] *C[j,k]* L2[k,l,c]*G[c,MM]*U[MM,M] + L2[j,i,a]*G[a,KK]*U[KK,K] *C[j,k]* L31[k,l,b,c]*G[c,MM]*U[MM,M]*G[b,LL]*V[LL,L])*G[l,JJ]*V[JJ,J];

        # fourth order tensor K4
        K4n[I,J,K,L] := 0.5*V[II,I]*G[i,II]*L1[j,i,k]*G[k,JJ]*V[JJ,J] *C[j,l]* L1[l,m,n]*G[n,KK]*V[KK,K]*G[m,LL]*V[LL,L];
        K5d[I,J,K,L,M] := V[II,I]*G[i,II]*(L1[j,i,a]*G[a,KK]*V[KK,K]*C[j,k]*L31[k,l,b,c]*G[c,LL]*U[LL,L]*G[b,MM]*V[MM,M] + L31[j,i,a,b]*G[b,LL]*U[LL,L]*G[a,KK]*V[KK,K]*C[j,k]*L1[k,l,c]*G[c,MM]*V[MM,M])*G[l,JJ]*V[JJ,J];
        K6dd[I,J,K,L,M,N] := 2.0*V[II,I]* G[i,II]*L31[j,i,a,b]*G[b,KK]*U[KK,K]*G[a,LL]*V[LL,L]* C[j,k]* L31[k,l,g,d]*G[d,MM]*U[MM,M]*G[g,NN]*V[NN,N]*G[l,JJ]*V[JJ,J];
    end
    K2n, K3d, K4dd, K3n, K4d, K5dd, K4n, K5d, K6dd
end
# STANDARD FORMULATION (no parametrization)
# Reduced Tensors at element-level (3d, evaluated at one gauss point)
function element_tensor_standard(G,C,H,L1,V)
    GHC = G'*H'*C
    CHG = C*H*G;
    @tensoropt begin
        K2n[I,J] := V[II,I]*G[i,II]*H[j,i]*C[j,k]*H[k,l]*G[l,JJ]*V[JJ,J];
        K3n[I,J,K] := 0.5*V[II,I]*GHC[II,i]*L1[i,j,l]*G[l,KK]*V[KK,K]*G[j,JJ]*V[JJ,J] + V[II,I]*G[i,II]*L1[j,i,k]*G[k,KK]*V[KK,K]*CHG[j,JJ]*V[JJ,J];
        K4n[I,J,K,L] := 0.5*V[II,I]*G[i,II]*L1[j,i,k]*G[k,JJ]*V[JJ,J] *C[j,l]* L1[l,m,n]*G[n,KK]*V[KK,K]*G[m,LL]*V[LL,L];
    end
    K2n, K3n, K4n
end


# STIFFNESS MATRIX and DERIVATIVES _____________________________________________
function stiffness_matrix_derivative(elements, nodes, connectivity, C, V, XGauss, WGauss)

    # useful dimensions
    nD = size(nodes,2)              # dimension of the problem
    nel = size(elements,1)          # number of elements
    nel_dofs = size(connectivity,2) # number of dofs per element
    nw = length(WGauss)             # number of integration points
    ndofs = prod(size(nodes))       # number of DOFs (total, free+constrained)

    G_FUN = Gselector(nD,nel_dofs)

    if nD==3
        L1, L2, L31, H = constantMatrices3D()
        A1, A2, A3 = Amatrices3D()
    elseif nD==2
        L1, L2, L31, H = constantMatrices2D()
        A1, A2, A3 = Amatrices2D()
    end

    # create indexes for sparse assembly
    I = reshape(kron(connectivity,ones(nel_dofs,1))',nel_dofs^2*nel,1)
    J = reshape(kron(connectivity,ones(1,nel_dofs))',nel_dofs^2*nel,1)
    I = I[:]
    J = J[:]
    dKdη_collection = zeros(nel_dofs^2, nel)

    # cycle over all the elements
    for e in 1:nel
        el_nodes = elements[e,:]        # IDs of the element nodes
        el_dofs  = connectivity[e,:]    # IDs of the element dofs
        xyz = nodes[el_nodes,:]         # element's coordinates x, y, z
        Ve = V[el_dofs]                 # Projection Basis V (only element's dofs)

        # cycle over Gauss quadrature points (integration over volume)
        for i in 1:nw
            ISO = XGauss[:,i]           # Gauss point coordinates
            w = WGauss[i]               # Gauss weight
            G, detJ = G_FUN(ISO,xyz)    # shape function derivative matrix G and Jacobian
            # compute the element-level dKdηE at the Gauss points
            CwJ = C .* (w * detJ);
            dKdηE = element_stiffness_matrix_derivative(G,CwJ,H,A1,L1,Ve)
            # assembly the element contributions in the global matrix
            dKdη_collection[:,e] .+= dKdηE[:]
        end
    end
    # sparse assembly
    dKdη_vectorized = reshape(dKdη_collection, nel_dofs^2*nel, 1)
    dKdη_vectorized = dKdη_vectorized[:]
    dKdη = sparse(I, J, dKdη_vectorized, ndofs, ndofs, +)
    # return all the tensors
    dKdη
end
function element_stiffness_matrix_derivative(G, C, H, A1, L1, Ve)
    th = G*Ve
    dKdη_E1 = G'*(H'*C*A1(th) + A1(th)'*C*H)*G
    CHGV = C*H*G*Ve;
    @tensoropt begin
        dKdη_E2[I,J] := G[ii,I] * ( L1[jj,ii,kk] * CHGV[jj] ) * G[kk,J]
    end
    dKdηE = dKdη_E1 + dKdη_E2
    dKdηE
end
function stiffness_matrix_sensitivity(elements, nodes, connectivity, C, U, XGauss, WGauss, formulation="N1")

    # useful dimensions
    nD = size(nodes,2)              # dimension of the problem
    nel = size(elements,1)          # number of elements
    nel_dofs = size(connectivity,2) # number of dofs per element
    nw = length(WGauss)             # number of integration points
    ndofs = prod(size(nodes))       # number of DOFs (total, free+constrained)

    G_FUN = Gselector(nD,nel_dofs)

    if nD==3
        L1, L2, L31, H = constantMatrices3D()
        A1, A2, A3 = Amatrices3D()
    elseif nD==2
        L1, L2, L31, H = constantMatrices2D()
        A1, A2, A3 = Amatrices2D()
    end

    if formulation == "N0"
        Ad = A1
    else
        Ad = A2 # for N1 and N1T
    end

    # create indexes for sparse assembly
    I = reshape(kron(connectivity,ones(nel_dofs,1))',nel_dofs^2*nel,1)
    J = reshape(kron(connectivity,ones(1,nel_dofs))',nel_dofs^2*nel,1)
    I = I[:]
    J = J[:]
    dKdξ_collection = zeros(nel_dofs^2, nel)

    # cycle over all the elements
    for e in 1:nel
        el_nodes = elements[e,:]        # IDs of the element nodes
        el_dofs  = connectivity[e,:]    # IDs of the element dofs
        xyz = nodes[el_nodes,:]         # element's coordinates x, y, z
        Ue = U[el_dofs]                 # Defect Basis U (only element's dofs)

        # cycle over Gauss quadrature points (integration over volume)
        for i in 1:nw
            ISO = XGauss[:,i]           # Gauss point coordinates
            w = WGauss[i]               # Gauss weight
            G, detJ = G_FUN(ISO,xyz)    # shape function derivative matrix G and Jacobian
            # compute the element-level dKdηE at the Gauss points
            CwJ = C .* (w * detJ);
            thd = G*Ue
            dKdξ_E = G'*(H'*CwJ*Ad(thd) + Ad(thd)'*CwJ*H)*G
            # assembly the element contributions in the global matrix
            dKdξ_collection[:,e] .+= dKdξ_E[:]
        end
    end
    # sparse assembly
    dKdξ_vectorized = reshape(dKdξ_collection, nel_dofs^2*nel, 1)
    dKdξ_vectorized = dKdξ_vectorized[:]
    dKdξ = sparse(I, J, dKdξ_vectorized, ndofs, ndofs, +)
    # return all the tensors
    dKdξ
end
function stiffness_matrix_linear(elements, nodes, connectivity, C, XGauss, WGauss)

    # useful dimensions
    nD = size(nodes,2)              # dimension of the problem
    nel = size(elements,1)          # number of elements
    nel_dofs = size(connectivity,2) # number of dofs per element
    nw = length(WGauss)             # number of integration points
    ndofs = prod(size(nodes))       # number of DOFs (total, free+constrained)

    G_FUN = Gselector(nD,nel_dofs)

    if nD==3
        L1, L2, L31, H = constantMatrices3D()
    elseif nD==2
        L1, L2, L31, H = constantMatrices2D()
    end

    # initialize tensors
    K0 = zeros(ndofs, ndofs)

    # cycle over all the elements
    for e in 1:nel
        el_nodes = elements[e,:]        # IDs of the element nodes
        el_dofs  = connectivity[e,:]    # IDs of the element dofs
        xyz = nodes[el_nodes,:]         # element's coordinates x, y, z

        # cycle over Gauss quadrature points (integration over volume)
        for i in 1:nw
            ISO = XGauss[:,i]           # Gauss point coordinates
            w = WGauss[i]               # Gauss weight
            G, detJ = G_FUN(ISO,xyz)    # shape function derivative matrix G and Jacobian
            # compute the element-level dKdηE at the Gauss points
            CwJ = C .* (w * detJ);
            K0E = G'*H'*CwJ*H*G
            # assembly the element contributions in the global matrix
            K0[el_dofs, el_dofs] .+= K0E
        end
    end
    # return all the tensors
    K0
end

# COMMONS & elements ___________________________________________________________
# all constant matrices are defined here (H,L₁,L₂,L₃)
function constantMatrices3D()
    # A = L��?⋅ θ
    L1 = zeros(6,9,9)
    L1[1,1,1]=1; L1[4,2,1]=1; L1[5,3,1]=1;
    L1[4,1,2]=1; L1[2,2,2]=1; L1[6,3,2]=1;
    L1[5,1,3]=1; L1[6,2,3]=1; L1[3,3,3]=1;
    L1[1,4,4]=1; L1[4,5,4]=1; L1[5,6,4]=1;
    L1[4,4,5]=1; L1[2,5,5]=1; L1[6,6,5]=1;
    L1[5,4,6]=1; L1[6,5,6]=1; L1[3,6,6]=1;
    L1[1,7,7]=1; L1[4,8,7]=1; L1[5,9,7]=1;
    L1[4,7,8]=1; L1[2,8,8]=1; L1[6,9,8]=1;
    L1[5,7,9]=1; L1[6,8,9]=1; L1[3,9,9]=1;

    # A₂ = L₂⋅ θᵈ
    L2 = zeros(6,9,9)
    L2[1,1,1]=1; L2[4,4,1]=1; L2[5,7,1]=1;
    L2[4,1,2]=1; L2[2,4,2]=1; L2[6,7,2]=1;
    L2[5,1,3]=1; L2[6,4,3]=1; L2[3,7,3]=1;
    L2[1,2,4]=1; L2[4,5,4]=1; L2[5,8,4]=1;
    L2[4,2,5]=1; L2[2,5,5]=1; L2[6,8,5]=1;
    L2[5,2,6]=1; L2[6,5,6]=1; L2[3,8,6]=1;
    L2[1,3,7]=1; L2[4,6,7]=1; L2[5,9,7]=1;
    L2[4,3,8]=1; L2[2,6,8]=1; L2[6,9,8]=1;
    L2[5,3,9]=1; L2[6,6,9]=1; L2[3,9,9]=1;
    L2 = -L2;

    # A₃A��? = (L₃��?⋅ θᵈ)⋅ θ
    L31 = zeros(6,9,9,9)
    L31[1,1,1,1] = 2; L31[4,2,1,1] = 1; L31[5,3,1,1] = 1; L31[4,1,2,1] = 1;
    L31[5,1,3,1] = 1; L31[1,4,4,1] = 2; L31[4,5,4,1] = 1; L31[5,6,4,1] = 1;
    L31[4,4,5,1] = 1; L31[5,4,6,1] = 1; L31[1,7,7,1] = 2; L31[4,8,7,1] = 1;
    L31[5,9,7,1] = 1; L31[4,7,8,1] = 1; L31[5,7,9,1] = 1; L31[4,1,1,2] = 2;
    L31[2,2,1,2] = 1; L31[6,3,1,2] = 1; L31[2,1,2,2] = 1; L31[6,1,3,2] = 1;
    L31[4,4,4,2] = 2; L31[2,5,4,2] = 1; L31[6,6,4,2] = 1; L31[2,4,5,2] = 1;
    L31[6,4,6,2] = 1; L31[4,7,7,2] = 2; L31[2,8,7,2] = 1; L31[6,9,7,2] = 1;
    L31[2,7,8,2] = 1; L31[6,7,9,2] = 1; L31[5,1,1,3] = 2; L31[6,2,1,3] = 1;
    L31[3,3,1,3] = 1; L31[6,1,2,3] = 1; L31[3,1,3,3] = 1; L31[5,4,4,3] = 2;
    L31[6,5,4,3] = 1; L31[3,6,4,3] = 1; L31[6,4,5,3] = 1; L31[3,4,6,3] = 1;
    L31[5,7,7,3] = 2; L31[6,8,7,3] = 1; L31[3,9,7,3] = 1; L31[6,7,8,3] = 1;
    L31[3,7,9,3] = 1; L31[1,2,1,4] = 1; L31[1,1,2,4] = 1; L31[4,2,2,4] = 2;
    L31[5,3,2,4] = 1; L31[5,2,3,4] = 1; L31[1,5,4,4] = 1; L31[1,4,5,4] = 1;
    L31[4,5,5,4] = 2; L31[5,6,5,4] = 1; L31[5,5,6,4] = 1; L31[1,8,7,4] = 1;
    L31[1,7,8,4] = 1; L31[4,8,8,4] = 2; L31[5,9,8,4] = 1; L31[5,8,9,4] = 1;
    L31[4,2,1,5] = 1; L31[4,1,2,5] = 1; L31[2,2,2,5] = 2; L31[6,3,2,5] = 1;
    L31[6,2,3,5] = 1; L31[4,5,4,5] = 1; L31[4,4,5,5] = 1; L31[2,5,5,5] = 2;
    L31[6,6,5,5] = 1; L31[6,5,6,5] = 1; L31[4,8,7,5] = 1; L31[4,7,8,5] = 1;
    L31[2,8,8,5] = 2; L31[6,9,8,5] = 1; L31[6,8,9,5] = 1; L31[5,2,1,6] = 1;
    L31[5,1,2,6] = 1; L31[6,2,2,6] = 2; L31[3,3,2,6] = 1; L31[3,2,3,6] = 1;
    L31[5,5,4,6] = 1; L31[5,4,5,6] = 1; L31[6,5,5,6] = 2; L31[3,6,5,6] = 1;
    L31[3,5,6,6] = 1; L31[5,8,7,6] = 1; L31[5,7,8,6] = 1; L31[6,8,8,6] = 2;
    L31[3,9,8,6] = 1; L31[3,8,9,6] = 1; L31[1,3,1,7] = 1; L31[4,3,2,7] = 1;
    L31[1,1,3,7] = 1; L31[4,2,3,7] = 1; L31[5,3,3,7] = 2; L31[1,6,4,7] = 1;
    L31[4,6,5,7] = 1; L31[1,4,6,7] = 1; L31[4,5,6,7] = 1; L31[5,6,6,7] = 2;
    L31[1,9,7,7] = 1; L31[4,9,8,7] = 1; L31[1,7,9,7] = 1; L31[4,8,9,7] = 1;
    L31[5,9,9,7] = 2; L31[4,3,1,8] = 1; L31[2,3,2,8] = 1; L31[4,1,3,8] = 1;
    L31[2,2,3,8] = 1; L31[6,3,3,8] = 2; L31[4,6,4,8] = 1; L31[2,6,5,8] = 1;
    L31[4,4,6,8] = 1; L31[2,5,6,8] = 1; L31[6,6,6,8] = 2; L31[4,9,7,8] = 1;
    L31[2,9,8,8] = 1; L31[4,7,9,8] = 1; L31[2,8,9,8] = 1; L31[6,9,9,8] = 2;
    L31[5,3,1,9] = 1; L31[6,3,2,9] = 1; L31[5,1,3,9] = 1; L31[6,2,3,9] = 1;
    L31[3,3,3,9] = 2; L31[5,6,4,9] = 1; L31[6,6,5,9] = 1; L31[5,4,6,9] = 1;
    L31[6,5,6,9] = 1; L31[3,6,6,9] = 2; L31[5,9,7,9] = 1; L31[6,9,8,9] = 1;
    L31[5,7,9,9] = 1; L31[6,8,9,9] = 1; L31[3,9,9,9] = 2;
    L31 ./= -2;

    # linear strain matrix H
    H = [1 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 1;
    0 1 0 1 0 0 0 0 0;
    0 0 1 0 0 0 1 0 0;
    0 0 0 0 0 1 0 1 0
    ];

    L1, L2, L31, H
end
function constantMatrices2D()
    # A = L��?⋅ θ
    L1 = zeros(3,4,4)
    L1[1,1,1] = 1; L1[3,2,1] = 1; L1[3,1,2] = 1; L1[2,2,2] = 1;
    L1[1,3,3] = 1; L1[3,4,3] = 1; L1[3,3,4] = 1; L1[2,4,4] = 1;

    # A₂ = L₂⋅ θᵈ
    L2 = zeros(3,4,4)
    L2[1,1,1] = 1; L2[3,3,1] = 1; L2[3,1,2] = 1; L2[2,3,2] = 1;
    L2[1,2,3] = 1; L2[3,4,3] = 1; L2[3,2,4] = 1; L2[2,4,4] = 1;
    L2 = -L2;

    # A₃A��? = (L₃��?⋅ θᵈ)⋅ θ
    L31 = zeros(3,4,4,4)
    L31[1,1,1,1] = 2; L31[3,2,1,1] = 1; L31[3,1,2,1] = 1; L31[1,3,3,1] = 2;
    L31[3,4,3,1] = 1; L31[3,3,4,1] = 1; L31[3,1,1,2] = 2; L31[2,2,1,2] = 1;
    L31[2,1,2,2] = 1; L31[3,3,3,2] = 2; L31[2,4,3,2] = 1; L31[2,3,4,2] = 1;
    L31[1,2,1,3] = 1; L31[1,1,2,3] = 1; L31[3,2,2,3] = 2; L31[1,4,3,3] = 1;
    L31[1,3,4,3] = 1; L31[3,4,4,3] = 2; L31[3,2,1,4] = 1; L31[3,1,2,4] = 1;
    L31[2,2,2,4] = 2; L31[3,4,3,4] = 1; L31[3,3,4,4] = 1; L31[2,4,4,4] = 2;
    L31 = -L31/2;

    # linear strain matrix H
    H = [1 0 0 0;
         0 0 0 1;
         0 1 1 0
    ];

    L1, L2, L31, H
end
# all function-matrices A are defined here (A₁,A₂,A₃)
function Amatrices3D()
    A1(th) =
        [th[1]     0     0 th[4]     0     0 th[7]     0     0;
             0 th[2]     0     0 th[5]     0     0 th[8]     0;
             0     0 th[3]     0     0 th[6]     0     0 th[9];
         th[2] th[1]     0 th[5] th[4]     0 th[8] th[7]     0;
         th[3]     0 th[1] th[6]     0 th[4] th[9]     0 th[7];
             0 th[3] th[2]     0 th[6] th[5]     0 th[9] th[8]]
    A2(th) = -1*
          [th[1]  th[4]  th[7]      0      0      0      0      0      0;
               0      0      0  th[2]  th[5]  th[8]      0      0      0;
               0      0      0      0      0      0  th[3]  th[6]  th[9];
           th[2]  th[5]  th[8]  th[1]  th[4]  th[7]      0      0      0;
           th[3]  th[6]  th[9]      0      0      0  th[1]  th[4]  th[7];
               0      0      0  th[3]  th[6]  th[9]  th[2]  th[5]  th[8]]
    A3(th) = -1/2*
        [2*th[1]       0       0       th[4]       th[7]           0;
               0 2*th[5]       0       th[2]           0       th[8];
               0	   0 2*th[9]           0       th[3]       th[6];
         2*th[2] 2*th[4]       0 th[1]+th[5]       th[8]       th[7];
         2*th[3]       0 2*th[7]       th[6] th[1]+th[9]       th[4];
             0   2*th[6] 2*th[8]       th[3]       th[2] th[5]+th[9]]
    A1, A2, A3
end
function Amatrices2D()
    A1(th) = [th[1]     0 th[3]     0;
                  0 th[2]     0 th[4];
              th[2] th[1] th[4] th[3]]
    A2(th) =-[th[1] th[3]     0     0;
                  0     0 th[2] th[4];
              th[2] th[4] th[1] th[3]]
    A3(th) = -1/2*
             [2*th[1]       0 th[3];
                    0 2*th[4] th[2];
              2*th[2] 2*th[3] th[1]+th[4]]
    A1, A2, A3
end
# shape functions selector
function Gselector(nDIM,ndofs)

    if nDIM == 3                # 3D continuum elements
        if ndofs == 24
            # 8-noded hexaedra
            G_FUN = G_HEX8
        elseif ndofs == 30
            # 10-noded tetrahedra
            G_FUN = G_TET10
        elseif ndofs == 45
            # 15-noded wedges
            G_FUN = G_WED15
        elseif ndofs == 60
            # 20-noded hexaedra
            G_FUN = G_HEX20
        end
    elseif nDIM == 2            # 2D continuum elements
        if ndofs == 8
            # 4-noded quadrilateral
            G_FUN = G_Q4
        elseif ndofs == 16
            # 8-noded quadrilateral
            G_FUN = G_Q8
        end
    end
    G_FUN
end
# shape functions derivatives for Q4 (4-noded quadrilateral)
function G_Q4(ISO,xy)
    r, s = ISO
    # shape function derivatives in natural coordinates
    dHn = 1/4*[s-1 1-s 1+s -1-s; r-1 -r-1 r+1 1-r];
    # jacobian
    J = dHn*xy;
    J1 = inv(J);
    detJ = det(J);
    # derivatives in physical coordinates
    dH = J1*dHn; # 3x4 matrix, [dNi_dx; dNi_dy] with i = 1...10
    # rearrange
    G = zeros(4,8)
    G[1:2,1:2:end] = dH;
    G[3:4,2:2:end] = dH;
    # return
    G, detJ
end
# shape functions derivatives for Q8 (8-noded quadrilateral)
function G_Q8(ISO,xy)
    g, h = ISO
    # shape function derivatives in natural coordinates
    dHn = 1/4*[-(2*g+h)*(h-1) -(h-1)*(2*g-h) (2*g+h)*(h+1) (h+1)*(2*g-h) 4*g*(h-1) 2-2*h^2 -4*g*(h+1) 2*h^2-2;
               -(g+2*h)*(g-1) -(g-2*h)*(g+1) (g+2*h)*(g+1) (g-2*h)*(g-1) 2*g^2-2 -4*h*(g+1) 2-2*g^2 4*h*(g-1)];
    # jacobian
    J = dHn*xy;
    J1 = inv(J);
    detJ = det(J);
    # derivatives in physical coordinates
    dH = J1*dHn; # 3x8 matrix, [dNi_dx; dNi_dy] with i = 1...10
    # rearrange
    G = zeros(4,16)
    G[1:2,1:2:end] = dH;
    G[3:4,2:2:end] = dH;
    # return
    G, detJ
end
# shape functions derivatives for TET10 (10-noded tetrahedra)
function G_TET10(ISO,xyz)
    g, h, r = ISO
    # shape function derivatives in natural coordinates
    dHn = zeros(3,10);
    dHn[1,1]=4*g+4*h+4*r-3; dHn[1,2]=4*g-1;         dHn[1,5]=4-4*h-4*r-8*g;
    dHn[1,6]=4*h;           dHn[1,7]=-4*h;          dHn[1,8]=-4*r;
    dHn[1,9]=4*r;           dHn[2,1]=4*g+4*h+4*r-3; dHn[2,3]=4*h-1;
    dHn[2,5]=-4*g;          dHn[2,6]=4*g;           dHn[2,7]=4-8*h-4*r-4*g;
    dHn[2,8]=-4*r;          dHn[2,10]=4*r;          dHn[3,1]=4*g+4*h+4*r-3;
    dHn[3,4]=4*r-1;         dHn[3,5]=-4*g;          dHn[3,7]=-4*h;
    dHn[3,8]=4-4*h-8*r-4*g; dHn[3,9]=4*g;           dHn[3,10]=4*h;
    # jacobian
    J = dHn*xyz;
    J1 = inv(J);
    detJ = det(J);
    # derivatives in physical coordinates
    dH = J1*dHn; # 3x10 matrix, [dNi_dx; dNi_dy; dNi_dz] with i = 1...10
    # rearrange
    G =  zeros(9,30)
    G[1:3,1:3:30] = dH;
    G[4:6,2:3:30] = dH;
    G[7:9,3:3:30] = dH;
    # return
    G, detJ
end
# shape functions derivatives for HEX8 (8-noded hexaedra)
function G_HEX8(ISO,xyz)
    g, h, r = ISO
    # shape function derivatives in natural coordinates
    dHn = [ -((h - 1)*(r - 1))/8  ((h - 1)*(r - 1))/8  -((h + 1)*(r - 1))/8  ((h + 1)*(r - 1))/8  ((h - 1)*(r + 1))/8  -((h - 1)*(r + 1))/8  ((h + 1)*(r + 1))/8  -((h + 1)*(r + 1))/8;
            -((g - 1)*(r - 1))/8  ((g + 1)*(r - 1))/8  -((g + 1)*(r - 1))/8  ((g - 1)*(r - 1))/8  ((g - 1)*(r + 1))/8  -((g + 1)*(r + 1))/8  ((g + 1)*(r + 1))/8  -((g - 1)*(r + 1))/8;
            -((g - 1)*(h - 1))/8  ((g + 1)*(h - 1))/8  -((g + 1)*(h + 1))/8  ((g - 1)*(h + 1))/8  ((g - 1)*(h - 1))/8  -((g + 1)*(h - 1))/8  ((g + 1)*(h + 1))/8  -((g - 1)*(h + 1))/8
    ];


    # jacobian
    J = dHn*xyz;
    J1 = inv(J);
    detJ = det(J);
    # derivatives in physical coordinates
    dH = J1*dHn; # 3x10 matrix, [dNi_dx; dNi_dy; dNi_dz] with i = 1...10
    # rearrange
    G =  zeros(9,24)
    G[1:3,1:3:24] = dH;
    G[4:6,2:3:24] = dH;
    G[7:9,3:3:24] = dH;
    # return
    G, detJ
end
# shape functions derivatives for HEX20 (20-noded hexaedra)
function G_HEX20(ISO,xyz)
    g, h, r = ISO
    # shape function derivatives in natural coordinates
    dHn=[((h-1)*(r-1)*(2*g+h+r+1))/8 -((h-1)*(r-1)*(h-2*g+r+1))/8 -((h+1)*(r-1)*(2*g+h-r-1))/8 -((h+1)*(r-1)*(2*g-h+r+1))/8 -((h-1)*(r+1)*(2*g+h-r+1))/8 -((h-1)*(r+1)*(2*g-h+r-1))/8 ((h+1)*(r+1)*(2*g+h+r-1))/8 ((h+1)*(r+1)*(2*g-h-r+1))/8 -(g*(h-1)*(r-1))/2 ((h^2-1)*(r-1))/4 (g*(h+1)*(r-1))/2 -((h^2-1)*(r-1))/4 (g*(h-1)*(r+1))/2 -((h^2-1)*(r+1))/4 -(g*(h+1)*(r+1))/2 ((h^2-1)*(r+1))/4 -((r^2-1)*(h-1))/4 ((r^2-1)*(h-1))/4 -((r^2-1)*(h+1))/4 ((r^2-1)*(h+1))/4;
         ((g-1)*(r-1)*(g+2*h+r+1))/8 -((g+1)*(r-1)*(2*h-g+r+1))/8 -((g+1)*(r-1)*(g+2*h-r-1))/8 -((g-1)*(r-1)*(g-2*h+r+1))/8 -((g-1)*(r+1)*(g+2*h-r+1))/8 -((g+1)*(r+1)*(g-2*h+r-1))/8 ((g+1)*(r+1)*(g+2*h+r-1))/8 ((g-1)*(r+1)*(g-2*h-r+1))/8 -((g^2-1)*(r-1))/4 (h*(g+1)*(r-1))/2 ((g^2-1)*(r-1))/4 -(h*(g-1)*(r-1))/2 ((g^2-1)*(r+1))/4 -(h*(g+1)*(r+1))/2 -((g^2-1)*(r+1))/4 (h*(g-1)*(r+1))/2 -((r^2-1)*(g-1))/4 ((r^2-1)*(g+1))/4 -((r^2-1)*(g+1))/4 ((r^2-1)*(g-1))/4;
         ((g-1)*(h-1)*(g+h+2*r+1))/8 -((g+1)*(h-1)*(h-g+2*r+1))/8 -((g+1)*(h+1)*(g+h-2*r-1))/8 -((g-1)*(h+1)*(g-h+2*r+1))/8 -((g-1)*(h-1)*(g+h-2*r+1))/8 -((g+1)*(h-1)*(g-h+2*r-1))/8 ((g+1)*(h+1)*(g+h+2*r-1))/8 ((g-1)*(h+1)*(g-h-2*r+1))/8 -((g^2-1)*(h-1))/4 ((h^2-1)*(g+1))/4 ((g^2-1)*(h+1))/4 -((h^2-1)*(g-1))/4 ((g^2-1)*(h-1))/4 -((h^2-1)*(g+1))/4 -((g^2-1)*(h+1))/4 ((h^2-1)*(g-1))/4 -(r*(g-1)*(h-1))/2 (r*(g+1)*(h-1))/2 -(r*(g+1)*(h+1))/2 (r*(g-1)*(h+1))/2
    ];

    # jacobian
    J = dHn*xyz;
    J1 = inv(J);
    detJ = det(J);
    # derivatives in physical coordinates
    dH = J1*dHn; # 3x10 matrix, [dNi_dx; dNi_dy; dNi_dz] with i = 1...10
    # rearrange
    G =  zeros(9,60)
    G[1:3,1:3:60] = dH;
    G[4:6,2:3:60] = dH;
    G[7:9,3:3:60] = dH;
    # return
    G, detJ
end
# shape functions derivatives for WED15 (15-noded wedges)
function G_WED15(ISO,xyz)
    g, h, r = ISO
    # shape function derivatives in natural coordinates
    dHn = zeros(3,15);
    dHn[1, 1]=1/2 - (r - 1)*(g + h - 1) - r^2/2 - ((r - 1)*(2*g + 2*h - 1))/2;
    dHn[1, 2]=r^2/2 - g*(r - 1) - ((2*g - 1)*(r - 1))/2 - 1/2;
    dHn[1, 4]=((r + 1)*(2*g + 2*h - 1))/2 + (r + 1)*(g + h - 1) - r^2/2 + 1/2;
    dHn[1, 5]=((2*g - 1)*(r + 1))/2 + g*(r + 1) + r^2/2 - 1/2;
    dHn[1, 7]=(r - 1)*(2*g + 2*h - 2) + 2*g*(r - 1);
    dHn[1, 8]=-2*h*(r - 1);
    dHn[1, 9]=2*h*(r - 1);
    dHn[1,10]=- (r + 1)*(2*g + 2*h - 2) - 2*g*(r + 1);
    dHn[1,11]=2*h*(r + 1);
    dHn[1,12]=-2*h*(r + 1);
    dHn[1,13]=r^2 - 1;
    dHn[1,14]=1 - r^2;
    dHn[2, 1]=1/2 - (r - 1)*(g + h - 1) - r^2/2 - ((r - 1)*(2*g + 2*h - 1))/2;
    dHn[2, 3]=r^2/2 - h*(r - 1) - ((2*h - 1)*(r - 1))/2 - 1/2;
    dHn[2, 4]=((r + 1)*(2*g + 2*h - 1))/2 + (r + 1)*(g + h - 1) - r^2/2 + 1/2;
    dHn[2, 6]=((2*h - 1)*(r + 1))/2 + h*(r + 1) + r^2/2 - 1/2;
    dHn[2, 7]=2*g*(r - 1);
    dHn[2, 8]=-2*g*(r - 1);
    dHn[2, 9]=2*(r - 1)*(g + h - 1) + 2*h*(r - 1);
    dHn[2,10]=-2*g*(r + 1);
    dHn[2,11]=2*g*(r + 1);
    dHn[2,12]=- 2*(r + 1)*(g + h - 1) - 2*h*(r + 1);
    dHn[2,13]=r^2 - 1;
    dHn[2,15]=1 - r^2;
    dHn[3, 1]=- ((g + h - 1)*(2*g + 2*h - 1))/2 - r*(g + h - 1);
    dHn[3, 2]=g*r - (g*(2*g - 1))/2;
    dHn[3, 3]=h*r - (h*(2*h - 1))/2;
    dHn[3, 4]=((g + h - 1)*(2*g + 2*h - 1))/2 - r*(g + h - 1);
    dHn[3, 5]=g*r + (g*(2*g - 1))/2;
    dHn[3, 6]=h*r + (h*(2*h - 1))/2;
    dHn[3, 7]=g*(2*g + 2*h - 2);
    dHn[3, 8]=-2*g*h;
    dHn[3, 9]=2*h*(g + h - 1);
    dHn[3,10]=-g*(2*g + 2*h - 2);
    dHn[3,11]=2*g*h;
    dHn[3,12]=-2*h*(g + h - 1);
    dHn[3,13]=2*r*(g + h - 1);
    dHn[3,14]=-2*g*r;
    dHn[3,15]=-2*h*r;

    # jacobian
    J = dHn*xyz;
    J1 = inv(J);
    detJ = det(J);
    # derivatives in physical coordinates
    dH = J1*dHn; # 3x10 matrix, [dNi_dx; dNi_dy; dNi_dz] with i = 1...10
    # rearrange
    G =  zeros(9,45)
    G[1:3,1:3:45] = dH;
    G[4:6,2:3:45] = dH;
    G[7:9,3:3:45] = dH;
    # return
    G, detJ
end

end
